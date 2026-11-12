// Publications on the home page, read live from the public ORCID record,
// so anything added to ORCID shows up here without editing the site.
(() => {
    const ORCID = '0009-0001-4199-8490';
    const API = `https://pub.orcid.org/v3.0/${ORCID}`;
    const ME = /wischnevsky/i;        // my name in the author lists
    const MAX_AUTHORS = 10;
    const list = document.getElementById('pubList');
    const count = document.getElementById('pubCount');
    const KINDS = { 'journal-article': 'Journal', preprint: 'Preprint', 'conference-paper': 'Conference', 'book-chapter': 'Chapter' };

    const getJSON = url => fetch(url, { headers: { Accept: 'application/json' } })
        .then(r => { if (!r.ok) throw new Error(r.status); return r.json(); });

    const value = x => (x && x.value) || '';

    function el(tag, cls, text) {
        const e = document.createElement(tag);
        if (cls) e.className = cls;
        if (text) e.textContent = text;
        return e;
    }

    function doi(work) {
        const ids = (work['external-ids'] || {})['external-id'] || [];
        const d = ids.find(i => i['external-id-type'] === 'doi');
        return d ? d['external-id-value'] : '';
    }

    // Journal name, or the arXiv number for a preprint
    function venue(work) {
        const journal = value(work['journal-title']);
        if (journal) return journal;
        const arxiv = (doi(work).match(/arxiv\.(.+)$/i) || value(work.url).match(/arxiv\.org\/abs\/(.+)$/) || [])[1];
        if (arxiv) return `arXiv:${arxiv}`;
        return work.type === 'preprint' ? 'Preprint' : '';
    }

    function link(work) {
        const url = value(work.url);
        if (url) return url;
        const d = doi(work);
        return d ? `https://doi.org/${d}` : '';
    }

    function date(work) {
        const d = work['publication-date'] || {};
        return [value(d.year), value(d.month) || '00', value(d.day) || '00'].join('-');
    }

    function authorsLine(work) {
        const names = ((work.contributors || {}).contributor || [])
            .map(c => value(c['credit-name'])).filter(Boolean);
        const p = el('p', 'pub-authors');
        const shown = names.slice(0, MAX_AUTHORS);
        shown.forEach((name, i) => {
            if (i) p.append(', ');
            p.append(ME.test(name) ? el('strong', '', name) : name);
        });
        if (names.length > shown.length) p.append(', et al.');
        return names.length ? p : null;
    }

    function render(works) {
        list.replaceChildren(...works.map(work => {
            const li = el('li', 'pub');
            const title = value((work.title || {}).title);
            const url = link(work);
            const t = el(url ? 'a' : 'span', 'pub-title', title);
            if (url) { t.href = url; t.target = '_blank'; t.rel = 'noopener'; }
            li.append(t);
            const authors = authorsLine(work);
            if (authors) li.append(authors);
            const meta = el('p', 'pub-meta');
            const kind = KINDS[work.type];
            if (kind) meta.append(el('span', `pub-badge ${work.type === 'journal-article' ? 'journal' : ''}`, kind));
            const year = value((work['publication-date'] || {}).year);
            meta.append(el('span', 'pub-venue', [venue(work), year].filter(Boolean).join(' · ')));
            li.append(meta);
            return li;
        }));
        count.textContent = works.length;
        count.hidden = false;
    }

    function fail() {
        list.replaceChildren(el('li', 'pub-note', "Couldn't load the list right now. See ORCID or Google Scholar above."));
    }

    // The summary has one entry per paper; the full records (with authors) come in one request
    getJSON(`${API}/works`)
        .then(data => {
            const codes = data.group.map(g => g['work-summary'][0]['put-code']).slice(0, 100);
            if (!codes.length) { list.replaceChildren(el('li', 'pub-note', 'None yet.')); return; }
            return getJSON(`${API}/works/${codes.join(',')}`).then(bulk => {
                const works = bulk.bulk.map(b => b.work).filter(Boolean);
                works.sort((a, b) => date(b).localeCompare(date(a)));
                render(works);
            });
        })
        .catch(fail);
})();
