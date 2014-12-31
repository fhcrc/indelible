
for i in /home/cwarth/src/matsen/indelible/help/example_files/*.shtml; do  
    echo $i
    f="$(basename $i)"
    f="${f%.*}"
    mkdir "$f"
    python ../html2text.py $i >$f/$f.txt
done
