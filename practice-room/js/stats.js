// Stats: progress by round, accuracy drawn on the staff, best and worst items, and mix-ups.
(() => {
  'use strict';
  const { $, CLEFS, draw } = App;
  const { H } = draw;

  // unit: what the time is per; item: what accuracy is tracked for; source: the drill
  // module that knows each item's label and notes (for drawing it on the staff);
  // levels: endless records are kept per difficulty; games: scored games instead of rounds
  const MODES = {
    clefs:        { name: 'Clefs',        unit: 'note',        item: 'note',     items: 'By note' },
    intervals:    { name: 'Intervals',    unit: 'interval',    item: 'interval', items: 'By interval', source: () => window.Intervals },
    chords:       { name: 'Chords',       unit: 'chord',       item: 'chord',    items: 'By chord',    source: () => window.Chords },
    melody:       { name: 'Melody',       unit: 'melody',                                                                          levels: true },
    progressions: { name: 'Progressions', unit: 'progression', item: 'chord',    items: 'By numeral',  source: () => window.Progressions, levels: true },
    intonation:   { name: 'Intonation',   unit: 'interval',    item: 'interval', items: 'By interval', source: () => window.Intonation,   levels: true },
    pitch:        { name: 'Pitch',        unit: 'note',        item: 'note',     items: 'By note',     source: () => window.Pitch },
    standards:    { name: 'Standards',    unit: 'song',                                                                            games: true },
  };
  const labelFor = (m, id) => {
    const src = MODES[m].source && MODES[m].source();
    return (src && src.label && src.label(id)) || id;
  };
  const MIN_TRIES = 3;                   // fewer tries: drawn faded, left out of best / worst
  const bin = acc => (acc >= 0.85 ? 'ok' : acc >= 0.6 ? 'mid' : 'bad');
  const pct = x => `${Math.round(x * 100)}%`;
  const secs = ms => `${(ms / 1000).toFixed(2)} s`;
  const sum = (arr, f) => arr.reduce((a, x) => a + f(x), 0);

  const screenEl = $('stats'), body = $('statsBody'), tip = $('statsTip');
  let mode = null, view = 'progress', clefShown = null;

  // ---------- Small DOM helpers (labels always go in as text) ----------
  function el(tag, cls, text) {
    const e = document.createElement(tag);
    if (cls) e.className = cls;
    if (text != null) e.textContent = text;
    return e;
  }
  function showTip(lines, x, y) {
    tip.replaceChildren(...lines.map((t, i) => el('div', i ? 'tip-sub' : 'tip-main', t)));
    tip.hidden = false;
    const r = tip.getBoundingClientRect();
    tip.style.left = `${Math.min(window.innerWidth - r.width - 8, Math.max(8, x + 14))}px`;
    tip.style.top = `${Math.max(8, y - r.height - 12)}px`;
  }
  const hideTip = () => { tip.hidden = true; };

  // ---------- Line chart (one series, by round) ----------
  function niceStep(max, count) {
    const raw = max / count, mag = 10 ** Math.floor(Math.log10(raw));
    return [1, 2, 2.5, 5, 10].map(m => m * mag).find(s => s >= raw);
  }

  // data: { pts: [{ y, tip: [lines] }], max?, fmt(v), empty }
  function drawLine(box) {
    const data = box.chartData;
    if (!data) return;
    const w = box.clientWidth, h = box.clientHeight;
    if (!w || !h) return;
    const { pts } = data;
    if (!pts.length) {
      box.replaceChildren(el('p', 'chart-empty', data.empty));
      return;
    }
    const m = { l: 40, r: 52, t: 10, b: 22 };
    const pw = Math.max(10, w - m.l - m.r), ph = Math.max(10, h - m.t - m.b);
    const top = data.max || (() => {
      const hi = Math.max(...pts.map(p => p.y)) * 1.1 || 1;
      const step = niceStep(hi, 3);
      return Math.ceil(hi / step) * step;
    })();
    const step = data.max ? data.max / 2 : niceStep(top, 3);
    const X = i => m.l + (pts.length === 1 ? pw / 2 : (i * pw) / (pts.length - 1));
    const Y = v => m.t + ph - (v / top) * ph;

    const parts = [];
    for (let v = 0; v <= top + 1e-9; v += step) {
      parts.push(`<line class="grid" x1="${m.l}" x2="${m.l + pw}" y1="${Y(v)}" y2="${Y(v)}"/>`);
      parts.push(`<text class="tick" x="${m.l - 8}" y="${Y(v) + 4}" text-anchor="end">${data.fmt(v, true)}</text>`);
    }
    const xTicks = [...new Set([0, Math.round((pts.length - 1) / 2), pts.length - 1])];
    for (const i of xTicks) parts.push(`<text class="tick" x="${X(i)}" y="${h - 6}" text-anchor="middle">${i + 1}</text>`);

    const line = pts.map((p, i) => `${i ? 'L' : 'M'}${X(i)} ${Y(p.y)}`).join('');
    if (pts.length > 1) {
      parts.push(`<path class="area" d="${line}L${X(pts.length - 1)} ${Y(0)}L${X(0)} ${Y(0)}Z"/>`);
      parts.push(`<path class="series" d="${line}"/>`);
    }
    const last = pts.length - 1;
    parts.push(`<circle class="dot" cx="${X(last)}" cy="${Y(pts[last].y)}" r="4"/>`);
    parts.push(`<text class="end-label" x="${X(last) + 9}" y="${Y(pts[last].y) + 4}">${data.fmt(pts[last].y)}</text>`);
    parts.push(`<g class="hover" visibility="hidden"><line class="crosshair" y1="${m.t}" y2="${m.t + ph}"/><circle class="dot" r="4"/></g>`);
    parts.push(`<rect class="hit" x="${m.l - 10}" y="0" width="${pw + 20}" height="${h}" fill="transparent"/>`);

    box.innerHTML = `<svg width="${w}" height="${h}" viewBox="0 0 ${w} ${h}" role="img">${parts.join('')}</svg>`;
    const svg = box.firstChild, hover = svg.querySelector('.hover');
    const move = e => {
      const r = svg.getBoundingClientRect();
      const i = pts.length === 1 ? 0 : Math.max(0, Math.min(last, Math.round(((e.clientX - r.left - m.l) / pw) * last)));
      hover.setAttribute('visibility', 'visible');
      hover.querySelector('line').setAttribute('x1', X(i));
      hover.querySelector('line').setAttribute('x2', X(i));
      hover.querySelector('circle').setAttribute('cx', X(i));
      hover.querySelector('circle').setAttribute('cy', Y(pts[i].y));
      showTip(pts[i].tip, e.clientX, e.clientY);
    };
    svg.addEventListener('pointermove', move);
    svg.addEventListener('pointerdown', move);
    svg.addEventListener('pointerleave', () => { hover.setAttribute('visibility', 'hidden'); hideTip(); });
  }

  // ---------- On the staff ----------
  // Clefs: a sideways histogram on the staff, one bar per line or space.
  function drawClefMap(box, clef, items) {
    const rows = new Map();
    for (const it of Object.values(items)) {
      if (!it.meta || it.meta.clef !== clef) continue;
      const r = rows.get(it.meta.p) || { n: 0, ok: 0, ms: 0, names: [] };
      r.n += it.n; r.ok += it.ok; r.ms += it.ms;
      r.names.push(`${it.meta.label}: ${it.ok} of ${it.n}`);
      rows.set(it.meta.p, r);
    }
    const c = CLEFS[clef];
    const yOf = p => (16.5 - p) * H, yBottom = yOf(0);
    const W = 330, vbH = yOf(-7.6);
    const noteX = 64, barX = 106, barL = 150, labelX = barX + barL + 8;
    const parts = [];
    for (let p = 0; p <= 8; p += 2) parts.push(draw.line('staff-line', 0, labelX - 4, yOf(p)));
    parts.push(draw.glyph('clef', c, 10, yBottom));
    for (const [p, r] of [...rows].sort((a, b) => b[0] - a[0])) {
      const acc = r.ok / r.n, y = yOf(p);
      const cls = `${bin(acc)}${r.n < MIN_TRIES ? ' faint' : ''}`;
      const nx = noteX + (p % 2 ? 24 : 0);   // neighbouring steps sit side by side, like a cluster
      parts.push(draw.ledgers(nx, p, yBottom));
      parts.push(draw.notePath(cls, nx, y));
      parts.push(`<rect class="map-track" x="${barX}" y="${y - 3}" width="${barL}" height="6" rx="3"/>`);
      parts.push(`<rect class="map-bar ${cls}" x="${barX}" y="${y - 3}" width="${Math.max(3, acc * barL)}" height="6" rx="3"/>`);
      parts.push(`<text class="map-label" x="${labelX}" y="${y + 2.6}">${App.noteName(c.bottomLine + p, 0)} · ${pct(acc)}</text>`);
      parts.push(`<rect class="hit" data-p="${p}" x="0" y="${y - H / 2}" width="${W}" height="${H}" fill="transparent"/>`);
    }
    box.innerHTML = `<svg class="map-svg" viewBox="0 0 ${W} ${vbH}" preserveAspectRatio="xMidYMid meet" role="img" aria-label="Accuracy for each note on the ${c.name.toLowerCase()} staff">${parts.join('')}</svg>`;
    box.firstChild.addEventListener('pointermove', e => {
      const p = e.target.dataset && e.target.dataset.p;
      if (p == null) return hideTip();
      const r = rows.get(Number(p));
      showTip([`${pct(r.ok / r.n)} correct · ${secs(r.ms / r.n)} per note`, ...r.names], e.clientX, e.clientY);
    });
    box.firstChild.addEventListener('pointerleave', hideTip);
  }

  // Intervals and chords: each drawn on middle C, with an accuracy column underneath.
  function drawColumnsMap(box, kind, items) {
    const mod = MODES[kind].source();
    const labelOf = id => labelFor(kind, id);
    const ids = mod.order.filter(id => items[id]);
    const colW = kind === 'chords' || kind === 'progressions' ? 58 : 46;
    const top = draw.trebleTop, bottom = top + 8 * H;
    const startX = 10 + CLEFS.treble.width + 14;
    const W = startX + ids.length * colW + 12;
    const barMax = 64, base = bottom + 100, vbY = 6, vbH = base + 22 - vbY;
    const parts = [];
    for (let p = 0; p <= 8; p += 2) parts.push(draw.line('staff-line', 0, W, bottom - p * H));
    parts.push(draw.glyph('clef', CLEFS.treble, 10, bottom));
    ids.forEach((id, i) => {
      const it = items[id], acc = it.ok / it.n;
      const cls = `${bin(acc)}${it.n < MIN_TRIES ? ' faint' : ''}`;
      const cx = startX + i * colW + colW / 2;
      parts.push(draw.columnParts(cx, mod.shape(id).map(n => ({ ...n, state: cls })), true));
      const bh = Math.max(3, acc * barMax), x = cx - 7, y = base - bh, r = Math.min(4, bh);
      parts.push(`<rect class="map-track" x="${x}" y="${base - barMax}" width="14" height="${barMax}" rx="4"/>`);
      parts.push(`<path class="map-bar ${cls}" d="M${x} ${base}V${y + r}Q${x} ${y} ${x + r} ${y}H${x + 14 - r}Q${x + 14} ${y} ${x + 14} ${y + r}V${base}Z"/>`);
      parts.push(`<text class="map-label" x="${cx}" y="${base - barMax - 5}" text-anchor="middle">${pct(acc)}</text>`);
      parts.push(`<text class="map-name" x="${cx}" y="${base + 15}" text-anchor="middle">${labelOf(id)}</text>`);
      parts.push(`<rect class="hit" data-id="${id}" x="${cx - colW / 2}" y="${vbY}" width="${colW}" height="${vbH}" fill="transparent"/>`);
    });
    // Fill the box; below a readable size, keep that size and scroll sideways
    const scale = Math.max(0.9, Math.min(box.clientHeight / vbH, box.clientWidth / W));
    box.innerHTML = `<svg class="map-svg scroll" width="${W * scale}" height="${vbH * scale}" viewBox="0 ${vbY} ${W} ${vbH}" role="img" aria-label="Accuracy for each ${MODES[kind].item}">${parts.join('')}</svg>`;
    box.firstChild.addEventListener('pointermove', e => {
      const id = e.target.dataset && e.target.dataset.id;
      if (!id) return hideTip();
      const it = items[id];
      showTip([`${labelOf(id)} · ${pct(it.ok / it.n)} correct`, `${it.ok} of ${it.n} · ${secs(it.ms / it.n)} each`], e.clientX, e.clientY);
    });
    box.firstChild.addEventListener('pointerleave', hideTip);
  }

  // ---------- Panels ----------
  function itemLabel(kind, id, it) {
    if (kind === 'clefs') return `${it.meta.label} · ${CLEFS[it.meta.clef].name.toLowerCase()}`;
    return labelFor(kind, id);
  }

  const tilesInto = (box, tiles) => box.replaceChildren(...tiles.map(([label, value]) => {
    const t = el('div', 'tile');
    t.append(el('span', 'tile-value', value), el('span', 'tile-label', label));
    return t;
  }));

  // Standard Guesser: high scores, then averages
  function renderGameTiles(s) {
    const g = s.games, n = g.length;
    const named = g.filter(x => x.songs);
    tilesInto($('statTiles'), [
      ['High score', n ? Math.max(...g.map(x => x.score)).toLocaleString('en-US') : '–'],
      ['Most songs', n ? String(Math.max(...g.map(x => x.songs))) : '–'],
      ['Games played', String(n)],
    ]);
    $('endlessTitle').textContent = 'Averages';
    tilesInto($('endlessTiles'), [
      ['Score', n ? Math.round(sum(g, x => x.score) / n).toLocaleString('en-US') : '–'],
      ['Songs / game', n ? (sum(g, x => x.songs) / n).toFixed(1) : '–'],
      ['Time per song', named.length ? secs(sum(named, x => x.avgMs * x.songs) / sum(named, x => x.songs)) : '–'],
    ]);
  }

  function renderGameCharts(s) {
    const empty = 'Play a game of Standard Guesser to start the chart.';
    $('timeTitle').textContent = 'Points, by game';
    $('accTitle').textContent = 'Songs named, by game';
    $('timeChart').chartData = {
      pts: s.games.map((g, i) => ({ y: g.score, tip: [`${g.score.toLocaleString('en-US')} points`, `Game ${i + 1} · ${g.songs} ${g.songs === 1 ? 'song' : 'songs'}`] })),
      fmt: (v, tick) => (tick && v >= 1000 ? `${+(v / 1000).toFixed(1)}k` : Math.round(v).toLocaleString('en-US')),
      empty,
    };
    $('accChart').chartData = {
      pts: s.games.map((g, i) => ({ y: g.songs, tip: [`${g.songs} ${g.songs === 1 ? 'song' : 'songs'} named`, `Game ${i + 1} · ${g.score.toLocaleString('en-US')} points`] })),
      fmt: v => String(Math.round(v)),
      empty,
    };
    drawLine($('timeChart'));
    drawLine($('accChart'));
  }

  function renderTiles(s) {
    if (MODES[mode].games) return renderGameTiles(s);
    const correct = sum(s.rounds, r => r.correct), total = sum(s.rounds, r => r.total);
    tilesInto($('statTiles'), [
      ['Rounds played', String(s.rounds.length)],
      ['Accuracy', total ? pct(correct / total) : '–'],
      [`Avg. time per ${MODES[mode].unit}`, s.rounds.length ? secs(sum(s.rounds, r => r.avgMs) / s.rounds.length) : '–'],
    ]);
    // Endless runs: the longest per difficulty, or for drills without one, the longest, how many and the average
    const runs = s.endless || {};
    if (MODES[mode].levels) {
      $('endlessTitle').textContent = 'Longest endless run';
      tilesInto($('endlessTiles'), [['easy', 'Easy'], ['medium', 'Medium'], ['hard', 'Hard']]
        .map(([k, name]) => [name, runs[k] ? String(runs[k].best) : '–']));
    } else {
      const r = runs.any;
      $('endlessTitle').textContent = 'Endless runs';
      tilesInto($('endlessTiles'), [
        ['Longest run', r ? String(r.best) : '–'],
        ['Runs played', String(r ? r.runs : 0)],
        ['Average run', r ? (r.total / r.runs).toFixed(1) : '–'],
      ]);
    }
  }

  function renderCharts(s) {
    if (MODES[mode].games) return renderGameCharts(s);
    const unit = MODES[mode].unit;
    $('accTitle').textContent = 'Accuracy, by round';
    const empty = `Finish a round of ${MODES[mode].name} to start the chart.`;
    $('timeTitle').textContent = `Avg. time per ${unit}, by round`;
    $('timeChart').chartData = {
      pts: s.rounds.map((r, i) => ({ y: r.avgMs / 1000, tip: [`${secs(r.avgMs)} per ${unit}`, `Round ${i + 1} · ${pct(r.correct / r.total)} correct`] })),
      fmt: (v, tick) => (tick ? `${+v.toFixed(2)}s` : `${v.toFixed(2)} s`),
      empty,
    };
    $('accChart').chartData = {
      pts: s.rounds.map((r, i) => ({ y: (100 * r.correct) / r.total, tip: [`${pct(r.correct / r.total)} correct`, `Round ${i + 1} · ${r.correct} of ${r.total}`] })),
      max: 100,
      fmt: v => `${Math.round(v)}%`,
      empty,
    };
    drawLine($('timeChart'));
    drawLine($('accChart'));
  }

  function renderMap(s) {
    const picker = $('clefPick');
    $('mapTitle').textContent = MODES[mode].items || '';
    const items = s.items;
    picker.hidden = mode !== 'clefs';
    if (mode === 'clefs') {
      const tries = k => sum(Object.values(items).filter(it => it.meta && it.meta.clef === k), it => it.n);
      if (!clefShown || !tries(clefShown)) clefShown = Object.keys(CLEFS).sort((a, b) => tries(b) - tries(a))[0];
      picker.replaceChildren(...Object.keys(CLEFS).map(k => {
        const b = el('button', null, CLEFS[k].name);
        b.type = 'button';
        b.setAttribute('role', 'radio');
        b.setAttribute('aria-checked', String(k === clefShown));
        b.disabled = !tries(k);
        b.addEventListener('click', () => { clefShown = k; renderMap(s); });
        return b;
      }));
    }
    renderBestWorst(s);
    drawMap(s);
  }

  function drawMap(s) {
    const box = $('mapBox'), items = s.items;
    if (!box.clientHeight) return;
    if (!Object.keys(items).length) {
      box.replaceChildren(el('p', 'chart-empty', `Play some ${MODES[mode].name} to see which ${MODES[mode].item}s you know best.`));
    } else if (mode === 'clefs') {
      drawClefMap(box, clefShown, items);
    } else {
      drawColumnsMap(box, mode, items);
    }
  }

  function renderBestWorst(s) {
    const list = Object.entries(s.items)
      .map(([id, it]) => ({ id, label: itemLabel(mode, id, it), n: it.n, acc: it.ok / it.n, ms: it.ms / it.n }))
      .filter(x => x.n >= MIN_TRIES);
    const box = $('bestWorst');
    if (list.length < 2) {
      box.replaceChildren(el('p', 'muted-note', `Best and worst ${MODES[mode].item}s appear after ${MIN_TRIES} tries each.`));
      return;
    }
    const best = [...list].sort((a, b) => b.acc - a.acc || a.ms - b.ms)[0];
    const worst = list.filter(x => x !== best).sort((a, b) => a.acc - b.acc || b.ms - a.ms).slice(0, 3);
    const row = (tag, x) => {
      const r = el('div', 'rank-row');
      const meter = el('span', 'meter');
      const fill = el('span', `meter-fill ${bin(x.acc)}`);
      fill.style.width = pct(x.acc);
      meter.append(fill);
      r.append(el('span', 'rank-tag', tag), el('span', 'rank-name', x.label), meter, el('span', 'rank-val', `${pct(x.acc)} · ${secs(x.ms)}`));
      return r;
    };
    box.replaceChildren(row('Best', best), ...worst.map((x, i) => row(i ? '' : 'Needs work', x)));
  }

  function renderMixups(s) {
    const box = $('mixups');
    const labelOf = id => labelFor(mode, id);
    const rows = Object.entries(s.mixups)
      .map(([k, n]) => { const [a, g] = k.split('>'); return { a, g, n }; })
      .sort((x, y) => y.n - x.n)
      .slice(0, 6);
    if (!rows.length) {
      box.replaceChildren(el('p', 'chart-empty', 'No mix-ups yet. When you pick the wrong answer, it shows up here.'));
      return;
    }
    const max = rows[0].n;
    box.replaceChildren(...rows.map(({ a, g, n }) => {
      const it = s.items[a];
      const misses = it ? it.n - it.ok : n;
      const r = el('div', 'mix-row');
      const name = el('div', 'mix-name');
      name.append(el('span', 'mix-chip', labelOf(a)), el('span', 'mix-as', 'as'), el('span', 'mix-chip', labelOf(g)));
      const bar = el('span', 'mix-bar');
      const fill = el('span', 'mix-fill');
      fill.style.width = `${(n / max) * 100}%`;
      bar.append(fill);
      r.title = `${labelOf(a)} answered as ${labelOf(g)} ${n} times: ${pct(Math.min(1, n / Math.max(1, misses)))} of your ${labelOf(a)} mistakes`;
      r.append(name, bar, el('span', 'mix-val', `${n}×`));
      return r;
    }), el('p', 'mix-note', 'Right answer first, then what you picked. Most frequent at the top.'));
  }

  // ---------- Screen ----------
  const viewsFor = m => (m === 'melody' || m === 'standards' ? ['progress'] : m === 'clefs' || m === 'intonation' ? ['progress', 'map'] : ['progress', 'map', 'mixups']);

  function render() {
    const s = App.readStats()[mode];
    for (const b of $('statsTabs').children) b.setAttribute('aria-selected', String(b.dataset.mode === mode));
    const views = viewsFor(mode);
    if (!views.includes(view)) view = 'progress';
    const viewBtns = $('statsViews');
    viewBtns.hidden = views.length < 2;
    for (const b of viewBtns.children) {
      b.hidden = !views.includes(b.dataset.view);
      b.setAttribute('aria-checked', String(b.dataset.view === view));
      if (b.dataset.view === 'map') b.textContent = MODES[mode].items || '';
    }
    body.dataset.mode = mode;
    body.dataset.view = view;
    hideTip();
    renderTiles(s);
    renderCharts(s);
    if (views.includes('map')) renderMap(s);
    if (views.includes('mixups')) renderMixups(s);
  }

  function open() {
    const all = App.readStats();
    if (!mode) {
      const last = m => Math.max(all[m].rounds.length ? all[m].rounds[all[m].rounds.length - 1].at : 0,
        ...Object.values(all[m].endless || {}).map(r => r.last || 0), ...(all[m].games || []).map(g => g.at));
      mode = Object.keys(MODES).sort((a, b) => last(b) - last(a))[0];
    }
    App.show('stats');
    render();
  }

  // ---------- Wiring ----------
  $('statsTabs').replaceChildren(...Object.entries(MODES).map(([id, m]) => {
    const b = el('button', null, m.name);
    b.type = 'button';
    b.setAttribute('role', 'tab');
    b.dataset.mode = id;
    b.addEventListener('click', () => { mode = id; render(); });
    return b;
  }));
  for (const b of $('statsViews').children) b.addEventListener('click', () => { view = b.dataset.view; render(); });

  let resetArmed = 0;
  $('statsReset').addEventListener('click', () => {
    const btn = $('statsReset');
    if (Date.now() - resetArmed < 4000) {
      App.resetStats();
      resetArmed = 0;
      btn.textContent = 'Reset stats';
      render();
      return;
    }
    resetArmed = Date.now();
    btn.textContent = 'Tap again to erase all stats';
    setTimeout(() => { if (Date.now() - resetArmed >= 4000) btn.textContent = 'Reset stats'; }, 4100);
  });

  // Redraw charts and maps whenever their boxes change size
  if (window.ResizeObserver) {
    let frame = 0;
    const dirty = new Set();
    const ro = new ResizeObserver(entries => {
      entries.forEach(e => dirty.add(e.target.id));
      cancelAnimationFrame(frame);
      frame = requestAnimationFrame(() => {
        if (App.current() !== 'stats' || !mode) return;
        for (const id of dirty) {
          if (id === 'mapBox') { if (viewsFor(mode).includes('map')) drawMap(App.readStats()[mode]); } else drawLine($(id));
        }
        dirty.clear();
      });
    });
    for (const id of ['timeChart', 'accChart', 'mapBox']) ro.observe($(id));
  }
  screenEl.addEventListener('scroll', hideTip, true);
  App.onRedraw(() => { if (App.current() === 'stats') render(); });

  window.Stats = { open };
})();
