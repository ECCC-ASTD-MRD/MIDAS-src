
outputFile='module_categories_new.txt'

# go into modules subdirectory
cd modules

rm $outputFile
rm module_category_list_raw.txt

for moduleFile in *_mod.f*
do
    for moduleName in `grep -i '^ *module ' $moduleFile |grep -iv "module *procedure"|cut -f2 -d' '`
    do
       category=`grep 'category=' $moduleFile |grep -i "$moduleName" | cut -f3 -d'='|sed "s/'//g"|cut -f3 -d'='|sed "s/)//g"`
       echo "FILE: $moduleFile"
       echo "  MODULE:   $moduleName"
       echo "  CATEGORY: $category" 
       echo $category >> module_category_list_raw.txt
    done
done

sort -u module_category_list_raw.txt > module_category_list.txt
rm module_category_list_raw.txt

IFS=$'\n'       # make newlines the only separator
for category in $(cat < module_category_list.txt); do
    echo $category >> $outputFile
    echo "================================" >> $outputFile
    for moduleFile in *_mod.f*
    do
	for moduleName in `grep 'category=' $moduleFile |grep -i "$category"|grep -o '\b\w*\_mod\b'`
	do
	    echo $moduleName >> $outputFile
	done
    done
    echo >> $outputFile
done

rm module_category_list.txt
