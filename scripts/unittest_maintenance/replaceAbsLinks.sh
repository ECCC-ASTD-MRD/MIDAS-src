#! /bin/bash

##-- Default values
ROOT=./
LINK_COMMON_ROOT=UnitTests
TARGET_COMMON_ROOT=UnitTests/midas
OUTPUT="::stdout::"
REPLACE=false

##-- Functions
function promptyn() {
  while true; do
    read -p "$1 [y/n]? " yn # Prompt the user and read input into 'yn'
    case $yn in
      [Yy]* ) return 0
      ;;
      [Nn]* ) return 1
      ;;
      * ) echo "Please answer yes or no."
      ;;
    esac
  done
}

##-- INTERFACE --
usage() {
echo "Find all symlinks in a file tree and identify those for wich the target is absolute."
echo "Under the assumption that the absolute target path share a branch leaf structure with"
echo "the root directory and that files are identical, identify how to replace the absolute"
echo "symlink target absolute path with a local relative path."
echo
echo "Usage: $0 [-d ROOT] [-s LINK_COMMON_ROOT] [-t TARGET_COMMON_ROOT] [-o OUTPUT] [-r -h]"
echo "  -d ROOT DIR           Specify the root directory (default: ./)"
echo "  -l LINK_COMMON_ROOT   Specify the common root ancestor for the symlink path (default: UnitTests)"
echo "  -t TARGET_COMMON_ROOT Specify the common root ancestor for it's target path (default: UnitTests/midas)"
echo "  -o OUTPUT             output file (default: stdout)"
echo "  -r                    proceed with actual replacement of absolute path with local relative path"
echo "  -h                    print this help"
}
## Parse command line options
while getopts "d:l:t:o:rh" opt; do
  case ${opt} in
    d )
      ROOT="$OPTARG"
      ;;
    l )
      LINK_COMMON_ROOT="$OPTARG"
      ;;
    t )
      TARGET_COMMON_ROOT="$OPTARG"
      ;;
    o )
      OUTPUT="$OPTARG"
      ;;
    r )
      REPLACE=true
      if promptyn "Do you really want to replace the symlinks?"; then
        continue
      else
        exit 1
      fi
      ;;
    h )
      usage; exit
      ;;
    * )
      usage; exit
      ;;
  esac
done

##-- MAIN --

declare -A LINKS_TARGET

ALL_LINKS=$(find ${ROOT} -type l -printf '%p:%l\n')

for duplet in ${ALL_LINKS}; do
  target=$(echo $duplet | cut -d':' -f 2)
  symlink=$(echo $duplet | cut -d':' -f 1)
  if [[ $target == /* ]]; then
    echo "$symlink:$target"

    # Extract directory part after 'UnitTests/'
    symlink_utests_part=$(echo "$symlink" | sed "s|.*/${LINK_COMMON_ROOT}/||")
    symlink_dir_utests=$(dirname "$symlink_utests_part")

    # Extract part after 'UnitTests/midas/'
    target_utests_part=$(echo "$target" | sed "s|.*/${TARGET_COMMON_ROOT}/||")
    target_dir_utests=$(dirname "$target_utests_part")
    target_filename=$(basename "$target_utests_part")

    # Split into arrays
    IFS='/' read -ra symlink_dir_parts <<< "$symlink_dir_utests"
    IFS='/' read -ra target_dir_parts <<< "$target_dir_utests"

    # Find common ancestor depth
    common_depth=0
    for ((i=0; i<${#symlink_dir_parts[@]} && i<${#target_dir_parts[@]}; i++)); do
      if [[ "${symlink_dir_parts[$i]}" == "${target_dir_parts[$i]}" ]]; then
        ((common_depth++))
      else
        break
      fi
    done

    # Calculate ups needed from symlink directory
    ups_needed=$((${#symlink_dir_parts[@]} - common_depth))
    ups=$(printf '../%.0s' $(seq 1 $ups_needed))

    # Extract relative part from target (after common ancestor)
    relative_path_parts=("${target_dir_parts[@]:$common_depth}")
    relative_path=$(printf '%s/' "${relative_path_parts[@]}")

    new_target="${ups}${relative_path}${target_filename}"

    LINKS_TARGET[$symlink]=$new_target
  fi
done

##-- OUTPUT AND REPLACEMENT --
for link in ${!LINKS_TARGET[@]}; do
  target=${LINKS_TARGET[$link]}
  cmd="rm -f $link && ln -s $target $link"
  if [ "${OUTPUT}" == "::stdout::" ]; then
    echo $cmd
  else
    echo $cmd >> ${OUTPUT}
  fi

  if ${REPLACE}; then
    eval $cmd
  fi
done
