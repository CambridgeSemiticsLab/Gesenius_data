from pathlib import Path

def compile_menu_items(dir, html_lines, level=0, base_dir=None):
    """Recursively compile analyses to match directory structures.."""
    for item in sorted(dir.glob('*')):
        if item.name == ".snakemake_timestamp": 
            continue
        indent = '&nbsp;'* 8 * level
        link = f'<a href="{item.relative_to(base_dir)}" target="_blank">{item.stem}</a>'
        html_lines.append(f'{indent}{link}')
        if item.is_dir():
            compile_menu_items(item, html_lines, level=level+1, base_dir=base_dir)

def make_html_menu(out_dir):
    """Build a menu to quickly access analyzed results."""
    out_dir = Path(out_dir)
    html_lines = []
    compile_menu_items(out_dir, html_lines, base_dir=out_dir)
    menu_html = '<html>{data}</html>'.format(data='\n<br>'.join(html_lines))
    out_dir.joinpath('menu.html').write_text(menu_html)

# construct the html menu
make_html_menu(snakemake.input.results)
