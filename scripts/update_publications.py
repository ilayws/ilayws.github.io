#!/usr/bin/env python3
"""Write publications.json from the public ORCID record (run daily by GitHub Actions)."""
import json
import re
import urllib.request
from pathlib import Path

ORCID = '0009-0001-4199-8490'
API = f'https://pub.orcid.org/v3.0/{ORCID}'
OUT = Path(__file__).resolve().parent.parent / 'publications.json'
KINDS = {'journal-article': 'Journal', 'preprint': 'Preprint', 'conference-paper': 'Conference', 'book-chapter': 'Chapter'}


def get(url):
    req = urllib.request.Request(url, headers={'Accept': 'application/json'})
    with urllib.request.urlopen(req, timeout=30) as r:
        return json.load(r)


def value(x):
    return (x or {}).get('value') or ''


def doi(work):
    ids = (work.get('external-ids') or {}).get('external-id') or []
    return next((i['external-id-value'] for i in ids if i['external-id-type'] == 'doi'), '')


def venue(work):
    journal = value(work.get('journal-title'))
    if journal:
        return journal
    m = re.search(r'arxiv\.(.+)$', doi(work), re.I) or re.search(r'arxiv\.org/abs/(.+)$', value(work.get('url')))
    if m:
        return f'arXiv:{m.group(1)}'
    return 'Preprint' if work.get('type') == 'preprint' else ''


def simplify(work):
    d = work.get('publication-date') or {}
    url = value(work.get('url')) or (f'https://doi.org/{doi(work)}' if doi(work) else '')
    return {
        'title': value((work.get('title') or {}).get('title')),
        'url': url,
        'authors': [value(c.get('credit-name')) for c in (work.get('contributors') or {}).get('contributor', []) if value(c.get('credit-name'))],
        'kind': KINDS.get(work.get('type'), ''),
        'venue': venue(work),
        'year': value(d.get('year')),
        'date': '-'.join([value(d.get('year')), value(d.get('month')) or '00', value(d.get('day')) or '00']),
    }


def main():
    codes = [g['work-summary'][0]['put-code'] for g in get(f'{API}/works')['group']]
    works = []
    for i in range(0, len(codes), 100):   # the bulk endpoint takes up to 100 at a time
        chunk = ','.join(str(c) for c in codes[i:i + 100])
        works += [simplify(b['work']) for b in get(f'{API}/works/{chunk}')['bulk'] if 'work' in b]
    works.sort(key=lambda w: w['date'], reverse=True)
    OUT.write_text(json.dumps(works, ensure_ascii=False, indent=2) + '\n', encoding='utf-8')
    print(f'{len(works)} publications -> {OUT.name}')


if __name__ == '__main__':
    main()
