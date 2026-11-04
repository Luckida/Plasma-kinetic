"""Check syntax and local links of the maintained GitHub documentation."""
from pathlib import Path
import ast
import json
import re

ROOT=Path(__file__).resolve().parent.parent
DOCS=['README.md','README_RU.md','CONTRIBUTING.md','CHANGELOG.md','NOTICE.md',
      'docs/PHYSICS_RU.md','docs/PHYSICS_CODE_MAP_RU.md','docs/REACTION_DATA_RU.md',
      'docs/ARCHITECTURE_RU.md','docs/TESTING_RU.md','docs/TROUBLESHOOTING_RU.md','docs/VALIDATION_2_8_2_RU.md','docs/DARK_HELIUM_2_8_3_RU.md']


def main():
    problems=[];files=list((ROOT/'src').glob('*.py'))+list((ROOT/'tools').glob('*.py'))+list((ROOT/'tests').glob('*.py'))+[ROOT/'run.py']
    for p in files:
        try:ast.parse(p.read_text(encoding='utf-8-sig'))
        except (SyntaxError,UnicodeError) as exc:problems.append(str(p.relative_to(ROOT))+': '+str(exc))
    for p in (ROOT/'configs').glob('*.json'):
        try:json.loads(p.read_text(encoding='utf-8-sig'))
        except (ValueError,UnicodeError) as exc:problems.append(str(p.relative_to(ROOT))+': '+str(exc))
    for name in DOCS:
        p=ROOT/name
        if not p.is_file():problems.append('Нет '+name);continue
        for link in re.findall(r'\]\(([^)]+)\)',p.read_text(encoding='utf-8')):
            if '://' in link or link.startswith('#'):continue
            target=link.split('#')[0]
            if target and not (p.parent/target).is_file():problems.append(f'{name}: нет цели {target}')
    if problems:raise SystemExit('\n'.join(problems))
    print(f'Python syntax: {len(files)} files; JSON and maintained local documentation links: OK.')


if __name__=='__main__':main()
