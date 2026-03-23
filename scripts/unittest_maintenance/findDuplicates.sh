#!/bin/bash

##-- Default values
ROOT=./
OUTPUT="::stdout::"
MAXSIZE_TO_PROCESS=0
N_THREADS=1
SEP=' '
REPLACE=false
DRYRUN=false



##-- INTERFACE --
usage() {
echo "Usage: $0 [-d ROOT] [-m MAXSIZE] [-n N_THREADS] [-o OUTPUT] [-r -t -h]"
echo "  -d ROOT DIR        Specify the root directory (default: ./)"
echo "  -m MAXSIZE         Maximum size of files to process (default: 0 - no maximum)"
echo "  -n N_THREADS       Number of threads to use (default: 1)"
echo "  -o OUTPUT          output file (default: stdout)"
echo "  -r                 replace duplicates with symlink"
echo "  -t                 test replacement only (dry run)"
echo "  -h                 print this help"
}

## Parse command line options
while getopts "d:m:n:o:rth" opt; do
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
    h )
      usage; exit
      ;;
    * )
      usage; exit
      ;;
  esac
done

##-- FUNCTIONS --

process_file() {
  local _file="$1"
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
  local _string_input="$1"
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

declare -A FILE_SIZE_ARRAY
declare -A FILE_HASH_ARRAY
declare -A HASH_LIST_ARRAY

export -f process_file
export ROOT
export MAXSIZE_TO_PROCESS

## Find files and process them in parallel, capturing output
results=$(find "$ROOT" -type f | xargs -P "$N_THREADS" -I {} bash -c 'process_file "$@"' _ {})

## Read results and populate associative arrays
while IFS=: read -r _file _size _md5sum; do
  FILE_SIZE_ARRAY["$_file"]="$_size"
  FILE_HASH_ARRAY["$_file"]="$_md5sum"
done <<< "$results"



## Triangular search
FILES=("${!FILE_HASH_ARRAY[@]}")
for ((i = 0; i < ${#FILES[@]}; i++)); do
  _fileA=${FILES[$i]}
  _hashA=${FILE_HASH_ARRAY[${_fileA}]}

  for ((j = i + 1; j < ${#FILES[@]}; j++)); do
    _fileB=${FILES[$j]}
    _hashB=${FILE_HASH_ARRAY[${_fileB}]}

    if [[ "$_hashA" == "$_hashB" ]]; then
      if [[ -z ${HASH_LIST_ARRAY[${_hashA}]} ]]; then
        HASH_LIST_ARRAY[${_hashA}]="${_fileA}${SEP}${_fileB}"
      else
        ## Split the existing entries into an array
        IFS=${SEP} read -r -a existing_files <<< "${HASH_LIST_ARRAY[${_hashA}]}"

        ## Append file if not already in the list
        if [[ ! " ${existing_files[@]} " =~ " ${_fileA} " ]]; then
          HASH_LIST_ARRAY[${_hashA}]+="${SEP}${_fileA}"
        fi
        if [[ ! " ${existing_files[@]} " =~ " ${_fileB} " ]]; then
          HASH_LIST_ARRAY[${_hashA}]+="${SEP}${_fileB}"
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

##-- OUTPUT AND REPLACEMENT --
if [ ! "${OUTPUT}" == "::stdout::" ]; then
  rm -f ${OUTPUT}
fi
for _hash in ${!HASH_LIST_ARRAY[@]}; do
  _string="> ${_hash}: ${HASH_LIST_ARRAY[${_hash}]}"
  if [ "${OUTPUT}" == "::stdout::" ]; then
    echo ${_string}
  else
    echo ${_string} >> ${OUTPUT}
  fi
  if ${REPLACE}; then
    replace_with_symlinks "${_string}"
  fi
done
