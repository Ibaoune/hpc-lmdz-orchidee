# This is a script in Bash.
# Author: Lionel GUEZ

# This script compares files "*.def" other than "*traceur.def" and
# "out.def" in two directories.

# This script uses GNU versions of the utilities cut, sort and uniq,
# and calls a Python script.

# See guide:
# http://lmdz.lmd.jussieu.fr/utilisateurs/outils/utilisation-de-lmdz#section-5

USAGE="usage: `basename $0` directory_1 directory_2"

if (($# != 2))
    then
    echo "$USAGE" >&2
    exit 1
fi

##set -x

type python3 >/dev/null 2>&1
if (($? == 0))
then
    my_python=python3
else
    my_python=python
fi

set -e

script_dir=`dirname $0`

# Let us make sure we know what the sort order is:
export LC_ALL=C

# Enable extended pattern matching features:
shopt -s extglob

for my_directory in $*
  do
  cd $my_directory
  rm -f all.def
  echo "def files in $my_directory:"
  ls !(*traceur|out).def

  # Concatenate and format the def files:
  cat !(*traceur|out).def | tr '\r' '\n' | tr -d " " \
      | grep -v -E '(^#|^!|INCLUDEDEF|^$)' \
      | tr '[:lower:]' '[:upper:]' | sed 's/=[NF]$/=FALSE/' \
      | sed 's/=[YT]$/=TRUE/' | sort >all.def

  # Check that no variable is repeated:
  cut --fields=1 --delimiter="=" all.def | uniq --repeated >plouf_$$
  if [[ -s plouf_$$ ]]
      then
      echo
      echo "Error: a variable is repeated in directory \"$my_directory\"."
      echo
      echo "Here is the list of repeated variables in directory \"$my_directory\":"
      echo
      cat plouf_$$
      rm plouf_$$
      rm all.def
      exit 1
  fi
  rm plouf_$$

  cd - >/dev/null
done


sort $1/all.def $2/all.def | uniq --repeated >common_lines

for my_directory in $*
  do
  sort $my_directory/all.def common_lines | uniq --unique \
      >$my_directory/uniq_lines
done

sort --field-separator="=" --key=1,1 --unique $2/uniq_lines $1/uniq_lines \
    | sort - $2/uniq_lines | uniq --unique >$1/uniq_var
sort --field-separator="=" --key=1,1 --unique $1/uniq_lines $2/uniq_lines \
    | sort - $1/uniq_lines | uniq --unique >$2/uniq_var

for my_directory in $*
  do
  cd $my_directory
  echo "Created file \"$my_directory/uniq_var\"."
  sort uniq_lines uniq_var | uniq --unique >common_var
  rm uniq_lines
  cd - >/dev/null
done

$my_python $script_dir/compare_real.py $*
sort common_lines --output=common_lines # was modified by the Python script
echo "Created files \"./common_lines\" and \"./uniq_val\"."

# Clean up:
for my_directory in $*
  do
    rm $my_directory/common_var $my_directory/all.def
done
