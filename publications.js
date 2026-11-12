// Publications on the home page. publications.json is refreshed from ORCID every day
// by .github/workflows/publications.yml, so new papers show up without editing the site.
(() => {
    const ME = /wischnevsky/i;        // my name in the author lists
    const MAX_AUTHORS = 10;
    const list = document.getElementById('pubList');
    const count = document.getElementById('pubCount');

    function el(tag, cls, text) {
        const e = document.createElement(tag);
        if (cls) e.className = cls;
        if (text) e.textContent = text;
        return e;
    }

    function authorsLine(names) {
        if (!names.length) return null;
        const p = el('p', 'pub-authors');
        const shown = names.slice(0, MAX_AUTHORS);
        shown.forEach((name, i) => {
            if (i) p.append(', ');
            p.append(ME.test(name) ? el('strong', '', name) : name);
        });
        if (names.length > shown.length) p.append(', et al.');
        return p;
    }

    function card(pub) {
        const li = el('li', 'pub');
        const t = el(pub.url ? 'a' : 'span', 'pub-title', pub.title);
        if (pub.url) { t.href = pub.url; t.target = '_blank'; t.rel = 'noopener'; }
        li.append(t);
        const authors = authorsLine(pub.authors || []);
        if (authors) li.append(authors);
        const meta = el('p', 'pub-meta');
        if (pub.kind) meta.append(el('span', pub.kind === 'Journal' ? 'pub-badge journal' : 'pub-badge', pub.kind));
        meta.append(el('span', 'pub-venue', [pub.venue, pub.year].filter(Boolean).join(' · ')));
        li.append(meta);
        return li;
    }

    function note(text) {
        list.replaceChildren(el('li', 'pub-note', text));
    }

    fetch('publications.json', { cache: 'no-cache' })
        .then(r => { if (!r.ok) throw new Error(r.status); return r.json(); })
        .then(pubs => {
            if (!pubs.length) return note('None yet.');
            list.replaceChildren(...pubs.map(card));
            count.textContent = pubs.length;
            count.hidden = false;
        })
        .catch(() => note("Couldn't load the list right now. See ORCID or Google Scholar above."));
})();
