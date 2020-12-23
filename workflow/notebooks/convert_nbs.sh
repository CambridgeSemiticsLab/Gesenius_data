OUTDIR="../notebooks_html"

for file in *.ipynb;
    do 
    jupyter nbconvert --to HTML --output-dir="$OUTDIR" $file;
    done
