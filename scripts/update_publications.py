#!/usr/bin/env python3
"""Write publications.json from ORCID, arXiv and Google Scholar (run by GitHub Actions).

Papers found in more than one place are merged by title. The journal version wins,
and its arXiv link is kept beside it. A source that can't be reached is skipped;
Google Scholar often refuses automated requests, so it only adds papers the others miss.
"""
import html
import json
import re
import sys
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from pathlib import Path

ORCID = '0009-0001-4199-8490'
SCHOLAR = '-En2ah0AAAAJ'
ARXIV_QUERY = 'au:Shlush_I OR au:Wischnevsky'     # "Shlush, I", "Shlush, I W", "Wischnevsky Shlush"
OUT = Path(__file__).resolve().parent.parent / 'publications.json'
UA = 'Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/141.0.0.0 Safari/537.36'
KINDS = {'journal-article': 'Journal', 'preprint': 'Preprint', 'conference-paper': 'Conference', 'book-chapter': 'Chapter'}
RANK = {'Journal': 0, 'Conference': 1, 'Chapter': 1, 'Preprint': 2, '': 3}


def fetch(url, accept='*/*'):
    req = urllib.request.Request(url, headers={'Accept': accept, 'User-Agent': UA, 'Accept-Language': 'en-US'})
    with urllib.request.urlopen(req, timeout=30) as r:
        return r.read().decode('utf-8', 'replace')


def is_me(name):
    """Shlush with first initial I (any spelling of the rest)."""
    return bool(re.search(r'shlush', name, re.I)) and name.strip()[:1].upper() == 'I'


def arxiv_id(text):
    m = re.search(r'(\d{4}\.\d{4,5})(v\d+)?', text or '')
    return m.group(1) if m else ''


def pub(title, url='', authors=(), kind='', venue='', date='', arxiv='', source=''):
    return {'title': ' '.join(title.split()), 'url': url, 'authors': list(authors), 'kind': kind,
            'venue': venue, 'year': date[:4], 'date': date, 'arxiv': arxiv, 'source': source}


# ---------- ORCID ----------
def orcid():
    api = f'https://pub.orcid.org/v3.0/{ORCID}'
    get = lambda url: json.loads(fetch(url, 'application/json'))
    value = lambda x: (x or {}).get('value') or ''
    codes = [g['work-summary'][0]['put-code'] for g in get(f'{api}/works')['group']]
    out = []
    for i in range(0, len(codes), 100):
        for b in get(f"{api}/works/{','.join(map(str, codes[i:i + 100]))}")['bulk']:
            w = b.get('work')
            if not w:
                continue
            ids = (w.get('external-ids') or {}).get('external-id') or []
            doi = next((x['external-id-value'] for x in ids if x['external-id-type'] == 'doi'), '')
            url = value(w.get('url')) or (f'https://doi.org/{doi}' if doi else '')
            ax = arxiv_id(doi) if 'arxiv' in doi.lower() else arxiv_id(url) if 'arxiv' in url else ''
            d = w.get('publication-date') or {}
            date = '-'.join(filter(None, [value(d.get('year')), value(d.get('month')), value(d.get('day'))]))
            kind = KINDS.get(w.get('type'), '')
            venue = value(w.get('journal-title')) or (f'arXiv:{ax}' if ax else '')
            authors = [value(c.get('credit-name')) for c in (w.get('contributors') or {}).get('contributor', [])]
            out.append(pub(value((w.get('title') or {}).get('title')), url, filter(None, authors), kind, venue, date, ax, 'orcid'))
    return out


# ---------- arXiv ----------
def arxiv():
    ns = {'a': 'http://www.w3.org/2005/Atom', 'x': 'http://arxiv.org/schemas/atom'}
    q = urllib.parse.urlencode({'search_query': ARXIV_QUERY, 'max_results': 200, 'sortBy': 'submittedDate'})
    root = ET.fromstring(fetch(f'https://export.arxiv.org/api/query?{q}'))
    out = []
    for e in root.findall('a:entry', ns):
        authors = [a.findtext('a:name', '', ns) for a in e.findall('a:author', ns)]
        if not any(is_me(a) for a in authors):
            continue
        ax = arxiv_id(e.findtext('a:id', '', ns))
        journal = (e.findtext('x:journal_ref', '', ns) or '').strip()
        doi = (e.findtext('x:doi', '', ns) or '').strip()
        url = f'https://doi.org/{doi}' if journal and doi else f'https://arxiv.org/abs/{ax}'
        out.append(pub(e.findtext('a:title', '', ns), url, authors, 'Journal' if journal else 'Preprint',
                       journal or f'arXiv:{ax}', e.findtext('a:published', '', ns)[:10], ax, 'arxiv'))
    return out


# ---------- Google Scholar ----------
def scholar():
    page = fetch(f'https://scholar.google.com/citations?user={SCHOLAR}&hl=en&cstart=0&pagesize=100&sortby=pubdate')
    rows = re.findall(r'<tr class="gsc_a_tr".*?</tr>', page, re.S)
    if not rows:
        raise RuntimeError('no rows (probably blocked)')
    text = lambda s: html.unescape(re.sub(r'<[^>]+>', '', s)).strip()
    out = []
    for row in rows:
        t = re.search(r'<a href="([^"]*)" class="gsc_a_at">(.*?)</a>', row, re.S)
        gray = [text(x) for x in re.findall(r'<div class="gs_gray">(.*?)</div>', row, re.S)]
        year = (re.search(r'gsc_a_h gsc_a_hc gs_ibl">(\d*)<', row) or [None, ''])[1]
        if not t:
            continue
        authors = [a.strip() for a in gray[0].split(',')] if gray else []
        venue = gray[1] if len(gray) > 1 else ''
        ax = arxiv_id(venue) if 'arxiv' in venue.lower() else ''
        url = f'https://arxiv.org/abs/{ax}' if ax else 'https://scholar.google.com' + html.unescape(t.group(1))
        venue = f'arXiv:{ax}' if ax else re.sub(r',\s*\d{4}$', '', venue)
        out.append(pub(text(t.group(2)), url, [a for a in authors if a != '...'], 'Preprint' if ax else '',
                       venue, year, ax, 'scholar'))
    return out


# ---------- Merge ----------
def key(title):
    return re.sub(r'[^a-z0-9]', '', title.lower())


def merge(found):
    groups = {}
    for p in found:
        groups.setdefault(key(p['title']), []).append(p)
    out = []
    for items in groups.values():
        # journal before preprint; ORCID and arXiv have full author lists, Scholar cuts them short
        items.sort(key=lambda p: (RANK.get(p['kind'], 3), p['source'] == 'scholar'))
        best = dict(items[0])
        best['authors'] = max((p['authors'] for p in items if p['source'] != 'scholar'), key=len, default=best['authors'])
        best['arxiv'] = next((p['arxiv'] for p in items if p['arxiv']), '')
        # the most precise date that agrees with the chosen version's
        best['date'] = max((p['date'] for p in items if p['date'].startswith(best['date'])), key=len)
        best['year'] = best['date'][:4]
        if best['kind'] == 'Journal' and best['arxiv']:
            best['also'] = f"https://arxiv.org/abs/{best['arxiv']}"
        if not best['kind']:
            best['kind'] = 'Preprint' if best['arxiv'] else ''
        best.pop('source')
        out.append(best)
    out.sort(key=lambda p: p['date'], reverse=True)
    return out


def main():
    found = []
    for name, source in (('ORCID', orcid), ('arXiv', arxiv), ('Google Scholar', scholar)):
        try:
            got = source()
            print(f'{name}: {len(got)}')
            found += got
        except Exception as e:   # one source down shouldn't wipe the list
            print(f'{name}: skipped ({e})', file=sys.stderr)
    if not found:
        raise SystemExit('No source could be reached; keeping the old list.')
    pubs = merge(found)
    OUT.write_text(json.dumps(pubs, ensure_ascii=False, indent=2) + '\n', encoding='utf-8')
    print(f'{len(pubs)} publications -> {OUT.name}')


if __name__ == '__main__':
    main()
