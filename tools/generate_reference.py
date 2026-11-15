#!/usr/bin/env python3
"""Regenerate source links and LaTeX/Markdown reference from maintained contracts.

Only text is read/written. Add a reviewed purpose entry to function_contracts.json
when adding a function; this tool does not invent documentation from a name.
"""
from pathlib import Path
import json,re

def prose(text):
    """Escape simple contract prose for LaTeX without treating it as commands."""
    return text.replace('&',r'\&').replace('%',r'\%').replace('_',r'\_').replace('#',r'\#').replace('^',r'\textasciicircum{}')

def main():
    """Validate all literal CLI options and refresh every maintained definition link."""
    r=Path(__file__).resolve().parents[1];inv=json.loads((r/'docs/function_contracts.json').read_text())
    pat=re.compile(r'(?m)^((?:static\s+|inline\s+|_Noreturn\s+)?(?:const\s+)?(?:unsigned\s+)?[A-Za-z_]\w*(?:\s*\*)?\s+(\*?[A-Za-z_]\w*)\s*\([^;{}]*?\))\s*\{')
    actual=[]
    for file in sorted((r/'src').glob('*.c')):
        for m in pat.finditer(file.read_text()):actual.append((str(file.relative_to(r)),m[2].lstrip('*')))
    assert set(actual)=={(x['file'],x['name']) for x in inv},'Review/add function contracts before generation'
    for x in inv:
        s=(r/x['file']).read_text();m=next(z for z in pat.finditer(s) if z[2].lstrip('*')==x['name']);x['line']=s[:m.start()].count('\n')+1;x['signature']=' '.join(m[1].split())
    (r/'docs/function_contracts.json').write_text(json.dumps(inv,indent=2)+'\n')
    md='# Runtime function reference\n\nEvery runtime C definition is listed below; see the [complete annotated-source atlas](../methodology/source-atlas.pdf).\n\n'
    tex=r'\section{Function-level implementation reference}'+'\nAll runtime C definitions are bound to reviewed responsibilities and actual source lines. Units and cache/ownership contracts accompany the source. The complete implementation, headers and tests are in the separate source atlas.\n'
    for file in sorted({x['file'] for x in inv}):
        group=[x for x in inv if x['file']==file];context=group[0]['module_contract']
        md+='## `'+file+'`\n\n'+context+'\n\n'
        tex+='\\subsection{\\texttt{'+file.replace('_',r'\_')+'}}\n'+prose(context)+'\n\\begin{description}[style=nextline,leftmargin=0pt,labelindent=0pt,labelwidth=0pt,labelsep=0pt]\n'
        for x in group:
            md+='* [`'+x['name']+'`](../'+file+'#L'+str(x['line'])+'): '+x['purpose']+'\n'
            tex+='\\item[\\texttt{'+x['name'].replace('_',r'\_')+'}] '+prose(x['purpose'])+' Source line '+str(x['line'])+'.\\par\n{\\small\\nolinkurl{'+x['signature']+'}}\n'
        md+='\n';tex+='\\end{description}\n'
    md+='## Inline ownership helpers\n\n`night_start(a)` is the global first-owned canonical offset. `night_owned(a)` is its owned-cell count without ghosts. Both are directly commented in `src/amr3_night_internal.h`.\n'
    tex+='\\subsection{Inline ownership helpers}\n\\code{night_start} returns the global molecular first-owned offset, and \\code{night_owned} its count excluding ghosts. Two canonical rows correspond to one internal row; k is contiguous. See the direct comments in \\path{src/amr3_night_internal.h}.\n'
    (r/'docs/FUNCTIONS.md').write_text(md);(r/'methodology/functions.tex').write_text(tex)
    items=json.loads((r/'docs/parameters.json').read_text());opts=set(re.findall(r'"(--[A-Za-z0-9-]+)"',(r/'src/main.c').read_text()));assert opts=={x['option'] for x in items}
    tex='\\section{Complete command-line reference}\nActual executable defaults are listed here; explicit teaching and qualification cases differ. Every option except help takes one value.\n\\begin{description}[style=nextline,leftmargin=0pt,labelindent=0pt,labelwidth=0pt,labelsep=0pt]\n'
    for x in items:tex+='\\item[\\texttt{'+prose(x['option'])+'}] \\textit{Type/units: '+prose(x['type_units'])+'; default: '+prose(x['default'])+'.} '+prose(x['description'])+'\n'
    tex+='\\end{description}\n';(r/'methodology/parameters.tex').write_text(tex)
    print('REFERENCE_GENERATION_PASS','functions='+str(len(inv)),'CLI_options='+str(len(items)))

if __name__=='__main__':main()
