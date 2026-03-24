#! /bin/bash

##-- Default values --
ROOT=./
MAXSIZE_TO_PROCESS=0
N_THREADS=1
SEP=' '
REPLACE=false
OUTPUT="::stdout::"
VERBOSE=false


##-- INTERFACE --
usage() {
echo "Usage: $0 [-d ROOT] [-m MAXSIZE] [-n N_THREADS] [-o OUTPUT] [-r -t -h -v]"
echo "  -d ROOT DIR        Specify the root directory (default: ./)"
echo "  -m MAXSIZE         Maximum size of files to process (default: 0 - no maximum)"
echo "  -n N_THREADS       Number of threads to use (default: 1)"
echo "  -o OUTPUT          output file (default: stdout)"
echo "  -r                 replace duplicates with symlink"
echo "  -t                 test replacement only (dry run)"
echo "  -v                 verbose outputs to stdout (output md5sum of all files)"
echo "  -h                 print this help"
}

## Parse command line options
while getopts "d:m:n:o:rthv" opt; do
  case ${opt} in
    d )
      ROOT="$OPTARG"
      ;;
    m )
      MAXSIZE_TO_PROCESS="$OPTARG"
      ;;
    n )
      N_THREADS="$OPTARG"
      ;;
    o )
      OUTPUT="$OPTARG"
      ;;
    r )
      REPLACE=true
      ;;
    t )
      DRYRUN=true
      ;;
    v )
      VERBOSE=true
      ;;
    h )
      usage; exit
      ;;
    * )
      usage; exit
      ;;
  esac
done

##-- FUNCTIONS --

function findTests () {
  ## find all unit test dir by looking for directories that only have leaf
  ## directories (unit test version directory)

  root=${1:?"ROOT!"}
  find $root -type d | while read dir; do
    subdir_count=$(find "$dir" -maxdepth 1 -type d | wc -l)
    [ "$subdir_count" -le 1 ] && continue

    leaf_count=$(find "$dir" -maxdepth 1 -type d -exec sh -c  '
      [ ! "$(find "$1" -maxdepth 1 -mindepth 1 -type d)" ] && echo 1
      ' sh {} \; | wc -l)
    [ "$((subdir_count - 1))" -eq "$leaf_count" ] && echo "$dir"
  done
}

process_file() {
  ## compute md5sum and size fo file
  local _file="${1:?"input file mandatory"}"
  local _md5sum

  local _size=$(stat --format="%s" "${_file}")

  if (( size <= MAXSIZE_TO_PROCESS )) || (( MAXSIZE_TO_PROCESS == 0 )); then
    local _md5sum
    _md5sum=$(md5sum "$_file" | cut -d' ' -f 1)
    echo "${_file}:${_size}:${_md5sum}"
  fi
}

relative_path() {
    # Arguments are the target file and the current working directory for the symlink
    set +x
    local _target="$1"
    local _current_dir="$2"

    # Get the absolute paths
    local _abs_target=$(realpath "$_target")
    local _abs_current_dir=$(realpath "$_current_dir")
    local _target_parts
    local _current_parts


    ## Convert absolute paths into arrays by splitting on /
    IFS='/' read -r -a _target_parts <<< ${_abs_target}
    IFS='/' read -r -a _current_parts <<< ${_abs_current_dir}
    ## the first element will be "" since the absolute path starts with "/"

    # Find the common base directory
    local _common_length=0
    for (( i = 0; i < ${#_target_parts[@]} && i < ${#_current_parts[@]}; i++ )); do
        if [[ "${_target_parts[i]}" == "${_current_parts[i]}" ]]; then
            common_length=$((i + 1))
        else
            break
        fi
    done

    # Build the relative path
    local rel_path=""
    # To construct the relative path, count how many directories to go back
    for (( i = common_length; i < ${#_current_parts[@]}; i++ )); do
        rel_path+="../"
    done
    for (( i = common_length; i < ${#_target_parts[@]}; i++ )); do
        rel_path+="${_target_parts[i]}"
        if [[ $i -lt $((${#_target_parts[@]} - 1)) ]]; then
            rel_path+="/"
        fi
    done

    echo "$rel_path"
}

replace_with_symlinks() {
  local _string_input="${1:?"string input mandatory"}"
  local _file_list
  local _kept_file
  local _string
  local _dup_dir
  local _dup_name
  local _dup_path
  local _link_rel_path
  local _duplicate
  local _path_array

  _file_list=$(echo ${_string_input} | cut -d':' -f 2 | xargs)
  IFS=' ' read -r -a _path_array <<< ${_file_list}
  if [ "${#_path_array[@]}" -lt 2 ]; then
    echo "Error: only one path, should not happen"
    exit 1
  fi
  ## we keep only the last file, other will be linked to it
  _kept_file="${_path_array[-1]}"
  unset '_path_array[-1]' ## poping it from the array

  ## replacing the other(s) with symlink to the kept one
  for _duplicate in ${_path_array[@]}; do
    _dup_dir=$(dirname ${_duplicate})
    _dup_name=$(basename ${_duplicate})
    _dup_path="${_dup_dir}/${_dup_name}"

    _link_rel_path=$(relative_path "${_kept_file}" "${_dup_dir}")

    _string="rm -f ${_duplicate} && ln -s ${_link_rel_path} ${_dup_path}"
    if [ "${OUTPUT}" == "::stdout::" ]; then
      echo ${_string}
    else
      echo ${_string} >> ${OUTPUT}
    fi
    if ! ${DRYRUN}; then
      rm -f ${_duplicate}
      ln -s ${_link_rel_path} ${_dup_path}
    fi
  done
}

##-- MAIN --

declare -A id_files_array
declare -A hash_in_test_array
declare -A file_size_array
declare -A file_hash_array

## Find all files and compute md5sum in parallel, capturing output
$VERBOSE && echo "Computing md5sum for all files in $ROOT"
export ROOT
export MAXSIZE_TO_PROCESS
export -f process_file
results=$(find "$ROOT" -type f | xargs -P "$N_THREADS" -I {} bash -c 'process_file "$@"' _ {})

## Read results and populate associative arrays
while IFS=: read -r _file _size _md5sum; do
  file_size_array["$_file"]="$_size"
  file_hash_array["$_file"]="$_md5sum"
done <<< "$results"


## find all unit tests
all_tests=$(findTests $ROOT)

for ut_test in $all_tests; do
  $VERBOSE && echo "Searching in $ut_test"

  ## Find all files within the different version of the test
  declare -A test_files
  n=0
  for file in $(find "$ut_test" -type f); do
    test_files[$n]=$file
    ((n++))
  done

  set +x ## listing clutter

  ## Triangular search within the test
  for ((i = 0; i < ${#test_files[@]}; i++)); do
    _fileA=${test_files[$i]}
    _hashA=${file_hash_array[${_fileA}]}

    $VERBOSE && echo "   ${_fileA/$ut_test/} : $_hashA"

    for ((j = i + 1; j < ${#test_files[@]}; j++)); do
      _fileB=${test_files[$j]}
      _hashB=${file_hash_array[${_fileB}]}

      if [[ "$_hashA" == "$_hashB" ]]; then
        if [[ -z ${id_files_array["$ut_test,${_hashA}"]} ]]; then
          hash_in_test_array[$ut_test]+=" $_hashA"
          id_files_array["$ut_test,${_hashA}"]="${_fileA}${SEP}${_fileB}"
        else
          ## Split the existing entries into an array
          IFS=${SEP} read -r -a existing_files <<< ${id_files_array["$ut_test,${_hashA}"]}

          ## Append file if not already in the list
          if [[ ! " ${existing_files[@]} " =~ " ${_fileA} " ]]; then
            id_files_array["$ut_test,${_hashA}"]+="${SEP}${_fileA}"
          fi
          if [[ ! " ${existing_files[@]} " =~ " ${_fileB} " ]]; then
            id_files_array["$ut_test,${_hashA}"]+="${SEP}${_fileB}"
          fi
        fi
      fi
    done

    ## Limit to N_THREADS
    if (( (i + 1) % N_THREADS == 0 )); then
      wait
    fi
  done

  ## Wait for any remaining background jobs
  wait
done


##-- OUTPUT AND REPLACEMENT --

if [ ! "${OUTPUT}" == "::stdout::" ]; then
  rm -f ${OUTPUT}
fi
for ut_test in $all_tests; do
  for _hash in ${hash_in_test_array[$ut_test]}; do
    _files=${id_files_array["$ut_test,${_hash}"]}
    _string="> ${_hash}: ${_files}"
    if [ "${OUTPUT}" == "::stdout::" ]; then
      echo ${_string//$ut_test/}
    else
      echo ${_string//$ut_test/} >> ${OUTPUT}
    fi
    if ${REPLACE}; then
      replace_with_symlinks "${_string}"
    fi
  done
done
