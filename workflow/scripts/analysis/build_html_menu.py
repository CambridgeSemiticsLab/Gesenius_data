import re
from pathlib import Path
import urllib

def sort_with_num(path):
    """Extract leading numbers in a file name for numerical sorting."""
    fname = path.name
    nums = re.match('^\d+', fname)
    if nums:
        return int(nums[0])
    else:
        return 0

get_indent = lambda level: '    ' * level 

def compile_menu_items(dir, html_lines, level=0, base_dir=None, outer=True):
    """Recursively compile analyses to match directory structures.."""

    # set up the outer list tag
    odent = get_indent(level)
    ul_class = ' class="tree"' if outer else ''
    html_lines.append(f'{odent}<ul{ul_class}>')

    # build the list items
    level += 1
    ident = get_indent(level)
    for item in sorted(dir.glob('*'), key=sort_with_num):

        # skip sys files
        if item.name.startswith('.') or item.name.endswith('.css'):
            continue

        # do directories recursively
        if item.is_dir():
            indent = get_indent(level)
            clean_name = re.match('\d+_(.*)', item.stem)
            if clean_name:
                name = clean_name.groups(0)[0]
            else:
                name = item.name
            html_lines.append(f'{ident}<li class="section">')
            html_lines.append(f'{ident}<input type="checkbox" id="{item}"/>')
            html_lines.append(f'{ident}<label for="{item}">{name}</label>')
            compile_menu_items(
                item, 
                html_lines, 
                level=level+1, 
                base_dir=base_dir, 
                outer=False
            )
            html_lines.append(f'{ident}</li>')

        # do files singly
        else:
            url = urllib.parse.quote(str(item.relative_to(base_dir)))
            link = f'<a href="{url}" target="_blank">{item.stem}</a>'
            li = f'{ident}<li>{link}</li>'
            html_lines.append(li)
        
    # finish the outer list
    html_lines.append(f'{odent}</ul>')

menu_doc = """\
<!DOCTYPE html>
<html>
    <head>
        <link rel="stylesheet" href="../../css/menu_style.css">
    </head>
    <body>
        <div class="menu_div">
{data}
        </div>
    </body>
</html>
"""

def make_html_menu(out_dir, level=0):
    """Build a menu to quickly access analyzed results."""
    out_dir = Path(out_dir)
    html_lines = []
    compile_menu_items(out_dir, html_lines, base_dir=out_dir, level=level)
    menu_html = menu_doc.format(data='\n'.join(html_lines))
    out_dir.joinpath('menu.html').write_text(menu_html)

# construct the html menu
make_html_menu(snakemake.params.results, level=3)
