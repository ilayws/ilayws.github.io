// Publications on the home page. publications.json is refreshed from ORCID, arXiv and
// Google Scholar by .github/workflows/publications.yml, so new papers show up on their own.
(() => {
    const ME = /shlush/i;             // my name in the author lists
    const MAX_AUTHORS = 10, SHOWN_BEFORE_ME = 6;
    const list = document.getElementById('pubList');
    const count = document.getElementById('pubCount');

    function el(tag, cls, text) {
        const e = document.createElement(tag);
        if (cls) e.className = cls;
        if (text) e.textContent = text;
        return e;
    }

    // Long lists are cut short, but always far enough to show my name
    function authorsLine(names) {
        if (!names.length) return null;
        const p = el('p', 'pub-authors');
        const add = name => p.append(ME.test(name) ? el('strong', '', name) : name);
        if (names.length <= MAX_AUTHORS) {
            names.forEach((name, i) => { if (i) p.append(', '); add(name); });
            return p;
        }
        const me = names.findIndex(n => ME.test(n));
        const head = me >= 0 && me < MAX_AUTHORS - 1 ? MAX_AUTHORS - 1 : SHOWN_BEFORE_ME;
        names.slice(0, head).forEach((name, i) => { if (i) p.append(', '); add(name); });
        if (me >= head) { p.append(', …, '); add(names[me]); }
        p.append(` et al. (${names.length} authors)`);
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
        if (pub.also) {
            const a = el('a', 'pub-also', 'arXiv');
            a.href = pub.also; a.target = '_blank'; a.rel = 'noopener';
            meta.append(a);
        }
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
