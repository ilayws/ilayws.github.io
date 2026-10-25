// Shared pieces: music helpers, staff drawing, piano audio, screens and results.
window.App = (() => {
  'use strict';

  // ---------- Music basics ----------
  // Notes are stored as diatonic numbers: octave * 7 + letter index (C4 = 28).
  const LETTERS = ['C', 'D', 'E', 'F', 'G', 'A', 'B'];
  const SEMIS = [0, 2, 4, 5, 7, 9, 11];
  const dia = (letter, octave) => octave * 7 + LETTERS.indexOf(letter);
  const diaName = d => LETTERS[d % 7] + Math.floor(d / 7);
  const diaMidi = d => 12 * (Math.floor(d / 7) + 1) + SEMIS[d % 7];

  // Staff position p = steps above the bottom line (0 = bottom line, 8 = top line).
  // top/bottom: where each clef's ink starts and ends, in the same steps.
  const CLEFS = {
    treble: { name: 'Treble', glyph: '\u{1D11E}', bottomLine: dia('E', 4), top: 10.8, bottom: -3.3 },
    bass:   { name: 'Bass',   glyph: '\u{1D122}', bottomLine: dia('G', 2), top: 8.05, bottom: 0.95 },
    alto:   { name: 'Alto',   glyph: '\u{1D121}', bottomLine: dia('F', 3), top: 8,    bottom: 0 },
  };
  // Accidental glyphs; top/bottom = ink extent in steps relative to the note.
  const ACCIDENTALS = {
    sharp:   { glyph: '♯', top: 2.8, bottom: -2.8 },
    flat:    { glyph: '♭', top: 3.5, bottom: -1.4 },
    natural: { glyph: '♮', top: 2.7, bottom: -2.7 },
    doubleSharp: { glyph: '\u{1D12A}', top: 1, bottom: -1 },
    doubleFlat:  { glyph: '\u{1D12B}', top: 3.5, bottom: -1.4 },
  };
  const ACC_NAMES = { 2: 'doubleSharp', 1: 'sharp', [-1]: 'flat', [-2]: 'doubleFlat' };
  const accFor = alter => ACC_NAMES[alter] || null;
  const ACC_SIGNS = { 2: '𝄪', 1: '♯', 0: '', [-1]: '♭', [-2]: '𝄫' };
  const noteName = (d, alter) => `${LETTERS[d % 7]}${ACC_SIGNS[alter] ?? ''}${Math.floor(d / 7)}`;
  const MUSIC_FONT = '"Noto Music", "Apple Symbols", "Segoe UI Symbol", serif';

  // ---------- Small helpers ----------
  const $ = id => document.getElementById(id);
  const randInt = (a, b) => a + Math.floor(Math.random() * (b - a + 1));
  const pick = arr => arr[Math.floor(Math.random() * arr.length)];
  const wait = ms => new Promise(r => setTimeout(r, ms));
  const reduceMotion = matchMedia('(prefers-reduced-motion: reduce)');

  function loadSettings(key, fallback, isValid) {
    try {
      const s = JSON.parse(localStorage.getItem(key));
      if (s && isValid(s)) return s;
    } catch (e) { /* storage unavailable */ }
    return fallback;
  }
  function saveSettings(key, value) {
    try { localStorage.setItem(key, JSON.stringify(value)); } catch (e) { /* ignore */ }
  }

  function makeTimers() {
    let list = [];
    return {
      later(fn, ms) { list.push(setTimeout(fn, ms)); },
      clear() { list.forEach(clearTimeout); list = []; },
    };
  }

  // ---------- Screens ----------
  let screen = 'hub';
  const GAME_SCREENS = new Set(['clefsGame', 'earGame', 'melodyGame', 'progGame']);
  const showHooks = [];
  function show(name) {
    if (GAME_SCREENS.has(screen) && name !== screen) stopAll();   // leaving a drill silences it
    const prev = screen;
    screen = name;
    for (const el of document.querySelectorAll('.screen')) el.hidden = el.id !== name;
    window.scrollTo(0, 0);
    showHooks.forEach(fn => fn(name, prev));
  }
  const current = () => screen;
  const onShow = fn => showHooks.push(fn);   // fn(newScreen, previousScreen)

  function setProgress(fillEl, textEl, done, total) {
    fillEl.closest('.progress').classList.remove('endless');
    fillEl.style.width = `${(Math.min(done, total) / total) * 100}%`;
    textEl.textContent = `${Math.min(done + 1, total)} / ${total}`;
  }

  // ---------- Endless runs ----------
  // Keep going until LIVES misses; the run's length is how many you got right.
  const LIVES = 3;
  const HEART = '<svg class="heart" viewBox="0 0 24 24" width="22" height="22" aria-hidden="true"><path d="M12 20.5s-7.6-4.7-9.6-9.4C1 7.7 3.3 4.5 6.7 4.5c2.2 0 3.9 1.3 5.3 3.2 1.4-1.9 3.1-3.2 5.3-3.2 3.4 0 5.7 3.2 4.3 6.6-2 4.7-9.6 9.4-9.6 9.4z"/></svg>';

  // The top bar during an endless run: hearts left instead of the progress bar, and the run
  // so far (or text, when given).
  function setRun(fillEl, textEl, lives, score, text) {
    const wrap = fillEl.closest('.progress');
    wrap.classList.add('endless');
    let hearts = wrap.querySelector('.hearts');
    if (!hearts) {
      hearts = document.createElement('span');
      hearts.className = 'hearts';
      hearts.setAttribute('role', 'img');
      hearts.innerHTML = HEART.repeat(LIVES);
      wrap.prepend(hearts);
    }
    [...hearts.children].forEach((h, i) => h.classList.toggle('lost', i >= lives));
    hearts.setAttribute('aria-label', `${lives} of ${LIVES} hearts left`);
    textEl.textContent = text ?? `Run ${score}`;
  }

  // The "per round" picker: the fixed counts, then Endless.
  // settings.count keeps the last number picked; settings.endless is true for Endless.
  // Returns a function that marks the current choice.
  function countPicker(el, counts, settings, onChange, idPrefix) {
    el.innerHTML = counts.map(n =>
      `<button type="button" role="radio" id="${idPrefix}-${n}" data-count="${n}" aria-checked="false">${n}</button>`).join('') +
      `<button type="button" role="radio" id="${idPrefix}-endless" class="endless-opt" data-count="endless" aria-checked="false" aria-label="Endless: keep going until you miss ${LIVES} times" title="Keep going until you miss ${LIVES} times"><span aria-hidden="true">∞</span><span class="endless-word" aria-hidden="true">Endless</span></button>`;
    el.addEventListener('click', e => {
      const b = e.target.closest('[data-count]');
      if (!b) return;
      settings.endless = b.dataset.count === 'endless';
      if (!settings.endless) settings.count = Number(b.dataset.count);
      onChange();
    });
    return () => {
      for (const b of el.children) {
        const on = settings.endless ? b.dataset.count === 'endless' : Number(b.dataset.count) === settings.count;
        b.setAttribute('aria-checked', String(on));
      }
    };
  }

  // ---------- Staff drawing (SVG) ----------
  const S = 16, H = S / 2;
  const yOf = p => (16 - p) * H;          // single staff: p = 16 sits at y = 0
  const FULL_RANGE = [15.4, -7.4];        // highest / lowest step a single staff shows
  const SLOT = 46, NOTE_RX = 12.5, NOTE_RY = 7.6, ACC_GAP = 4;

  // Fit each clef and accidental glyph to its spot on the staff by measuring
  // its ink, so it lines up whichever music font the browser ends up using.
  function measureGlyphs() {
    const ctx = document.createElement('canvas').getContext('2d');
    ctx.font = `100px ${MUSIC_FONT}`;
    for (const g of [...Object.values(CLEFS), ...Object.values(ACCIDENTALS)]) {
      const m = ctx.measureText(g.glyph);
      const inkH = m.actualBoundingBoxAscent + m.actualBoundingBoxDescent;
      if (!inkH) { Object.assign(g, { size: 40, rise: 0, dx: 0, width: 12 }); continue; }
      const k = ((g.top - g.bottom) * H) / inkH;
      g.size = 100 * k;
      g.rise = m.actualBoundingBoxAscent * k - g.top * H;   // baseline = anchor y + rise
      g.dx = m.actualBoundingBoxLeft * k;
      g.width = (m.actualBoundingBoxLeft + m.actualBoundingBoxRight) * k;
    }
  }

  function ellipsePath(cx, cy, a, b, deg) {
    const t = deg * Math.PI / 180, ux = a * Math.cos(t), uy = a * Math.sin(t);
    const x1 = cx - ux, y1 = cy - uy, x2 = cx + ux, y2 = cy + uy;
    return `M${x1} ${y1}A${a} ${b} ${deg} 1 0 ${x2} ${y2}A${a} ${b} ${deg} 1 0 ${x1} ${y1}Z`;
  }
  // Whole note: an oval with a tilted hole (even-odd fill).
  const wholeNote = (cx, cy) => ellipsePath(cx, cy, NOTE_RX, NOTE_RY, 0) + ellipsePath(cx, cy, 6, 3.9, -55);
  const notePath = (cls, cx, cy) => `<path class="note ${cls || ''}" fill-rule="evenodd" d="${wholeNote(cx, cy)}"/>`;

  const line = (cls, x1, x2, y) => `<line class="${cls}" x1="${x1}" x2="${x2}" y1="${y}" y2="${y}"/>`;
  const glyph = (cls, g, x, y) =>
    `<text class="${cls}" x="${x + g.dx}" y="${y + g.rise}" font-size="${g.size}" font-family='${MUSIC_FONT}'>${g.glyph}</text>`;

  // Ledger lines for a note at step p, on a staff whose bottom line is at y = yBottom.
  function ledgers(cx, p, yBottom) {
    const out = [];
    for (let q = -2; q >= p; q -= 2) out.push(line('ledger', cx - NOTE_RX - 5, cx + NOTE_RX + 5, yBottom - q * H));
    for (let q = 10; q <= p; q += 2) out.push(line('ledger', cx - NOTE_RX - 5, cx + NOTE_RX + 5, yBottom - q * H));
    return out.join('');
  }

  // Scale to the available width (within min/max), then shrink further if
  // the drawing would be taller than maxHeight, so the screen never scrolls.
  function sizeSvg(svg, x, y, w, h, opts) {
    svg.setAttribute('viewBox', `${x} ${y} ${w} ${h}`);
    let scale = opts.scale ?? Math.min(opts.max ?? 1.5, Math.max(opts.min ?? 1, (opts.fit || w) / w));
    if (opts.maxHeight) scale = Math.max(0.4, Math.min(scale, opts.maxHeight / h));
    svg.style.width = `${w * scale}px`;
    svg.style.height = `${h * scale}px`;
    return scale;
  }

  /**
   * Single staff with clef changes.
   * items: [{type:'clef', clef}, {type:'note', p, acc, index, state}, {type:'chord', ps}]
   * opts: fit (available px), min/max scale, barline, range [topStep, bottomStep]
   * Returns { scale, noteX } so callers can scroll to a note.
   */
  function renderStaff(svg, items, opts = {}) {
    const body = [];
    const noteX = [];
    const [rTop, rBottom] = opts.range || FULL_RANGE;
    const yBottom = yOf(0);
    let x = 10;

    for (const it of items) {
      if (it.type === 'clef') {
        const c = CLEFS[it.clef];
        x += 8;
        body.push(glyph('clef', c, x, yBottom));
        x += c.width + 14;
      } else if (it.type === 'chord') {
        const cx = x + SLOT / 2;
        for (const p of it.ps) body.push(ledgers(cx, p, yBottom) + notePath('', cx, yOf(p)));
        x += SLOT;
      } else {
        const acc = it.acc && ACCIDENTALS[it.acc];
        const accW = acc ? acc.width + ACC_GAP : 0;
        const cx = x + accW + SLOT / 2;
        const y = yOf(it.p);
        noteX[it.index] = cx;
        if (it.state === 'current') {
          body.push(`<rect class="cursor" x="${x + 3}" y="${yOf(rTop - 0.2)}" width="${accW + SLOT - 6}" height="${yOf(rBottom + 0.2) - yOf(rTop - 0.2)}" rx="9"/>`);
        }
        body.push(ledgers(cx, it.p, yBottom));
        if (acc) body.push(glyph(`acc ${it.state || ''}`, acc, cx - NOTE_RX - ACC_GAP - acc.width, y));
        body.push(notePath(it.state, cx, y));
        x += accW + SLOT;
      }
    }

    let end;
    const lines = [];
    if (opts.barline === false) {
      end = x + 16;
    } else {
      end = x + 18;
      lines.push(`<line class="barline" x1="${end - 9}" x2="${end - 9}" y1="${yOf(8)}" y2="${yOf(0)}"/>`);
      lines.push(`<rect class="barline-thick" x="${end - 5}" y="${yOf(8)}" width="4.5" height="${yOf(0) - yOf(8)}"/>`);
      end -= 0.5;
    }
    for (let p = 0; p <= 8; p += 2) lines.unshift(line('staff-line', 0, end, yOf(p)));

    svg.innerHTML = lines.join('') + body.join('');
    const vbY = yOf(rTop);
    const scale = sizeSvg(svg, 0, vbY, end, yOf(rBottom) - vbY, opts);
    return { scale, noteX };
  }

  // ---------- Grand staff (treble + bass) ----------
  const GS_TOP = 46;                       // y of the treble top line (room for two ledger lines above)
  const GS_T_BOTTOM = GS_TOP + 8 * H;      // treble bottom line
  const GS_B_TOP = GS_T_BOTTOM + 8 * H;    // bass top line
  const GS_B_BOTTOM = GS_B_TOP + 8 * H;    // bass bottom line
  const GS_LABEL_Y = GS_B_BOTTOM + 36;
  const GS_HEIGHT = GS_LABEL_Y + 14;
  const MIDDLE_C = dia('C', 4);

  function bracePath(x, top, bottom) {
    const h = bottom - top, m = (top + bottom) / 2, w = 11;
    return `M${x + w} ${top}` +
      `C${x + 0.5} ${top + h * 0.05} ${x + w * 0.9} ${m - h * 0.12} ${x} ${m}` +
      `C${x + w * 0.9} ${m + h * 0.12} ${x + 0.5} ${bottom - h * 0.05} ${x + w} ${bottom}` +
      `C${x + w * 0.8} ${bottom - h * 0.08} ${x + w * 2.1} ${m + h * 0.1} ${x} ${m}` +
      `C${x + w * 2.1} ${m - h * 0.1} ${x + w * 0.8} ${top + h * 0.08} ${x + w} ${top}Z`;
  }

  // Height available to a drawing inside a card that sits in a flex "stage":
  // the stage's height minus everything in the card that isn't the drawing.
  function roomFor(stage, card, drawingBox) {
    return stage.clientHeight - (card.offsetHeight - drawingBox.offsetHeight) - 2;
  }

  // Draw with render(maxHeight), then check the card really fits its stage;
  // if anything around it changed size, shrink the drawing by the overflow.
  function fitInStage(stage, card, drawingBox, render) {
    const layout = render(roomFor(stage, card, drawingBox));
    const over = card.offsetHeight - stage.clientHeight;
    return over > 0 ? render(drawingBox.offsetHeight - over - 2) : layout;
  }

  // Redraw whenever the stage's size changes (window resize, fonts loading, browser bars).
  function watchStage(stage, redrawFn) {
    if (!window.ResizeObserver) return;
    let frame = 0;
    new ResizeObserver(() => {
      cancelAnimationFrame(frame);
      frame = requestAnimationFrame(redrawFn);
    }).observe(stage);
  }

  // One column of notes (a chord, or a single note) centred on x = center.
  // Seconds sit side by side; accidentals step left when they would collide.
  function columnParts(center, colNotes, single) {
    const out = [];
    // Place each note: staff, step, x offset (seconds sit side by side)
    const notes = colNotes.map(n => {
      const treble = single || n.d >= MIDDLE_C;
      const yBottom = treble ? GS_T_BOTTOM : GS_B_BOTTOM;
      const p = n.d - (treble ? CLEFS.treble : CLEFS.bass).bottomLine;
      return { ...n, treble, yBottom, p, y: yBottom - p * H, off: 0 };
    }).sort((a, b) => a.d - b.d);
    for (let i = 1; i < notes.length; i++) {
      const lo = notes[i - 1], hi = notes[i];
      if (lo.treble === hi.treble && hi.p - lo.p <= 1 && lo.off === 0) hi.off = NOTE_RX * 2 - 1;
    }
    // Accidentals, top note first, stepping left when two would collide
    const placed = [];
    let left = -NOTE_RX;
    for (const n of [...notes].reverse()) {
      const acc = ACCIDENTALS['acc' in n ? n.acc : accFor(n.alter)];
      if (!acc) continue;
      let right = -NOTE_RX - ACC_GAP;
      for (const q of placed) if (Math.abs(q.y - n.y) < 5.5 * H) right = Math.min(right, q.x - 2);
      const ax = right - acc.width;
      placed.push({ y: n.y, x: ax, acc, n });
      left = Math.min(left, ax);
    }
    const rightEdge = Math.max(...notes.map(n => n.off)) + NOTE_RX;
    const shift = center - (left + rightEdge) / 2;

    for (const n of notes) {
      const cx = shift + n.off;
      out.push(ledgers(cx, n.p, n.yBottom));
      out.push(notePath(n.state, cx, n.y));
    }
    for (const a of placed) out.push(glyph(`acc ${a.n.state || ''}`, a.acc, shift + a.x, a.y));
    return out.join('');
  }

  /**
   * columns: array (length = opts.cols) of null or
   *   { notes: [{d, alter, acc?, state}], label?, labelState?, cursor? }
   * Notes from middle C up go on the treble staff, the rest on the bass staff.
   * col.sublabel adds a second label line (opts.sublabels makes room for it).
   * A note's acc overrides the accidental implied by its alter ('' = none).
   * opts.trebleOnly draws a single treble staff that takes every note.
   */
  function renderGrandStaff(svg, columns, opts = {}) {
    const cols = opts.cols || columns.length;
    const colW = opts.colWidth || 84;
    const single = !!opts.trebleOnly;
    const sysBottom = single ? GS_T_BOTTOM : GS_B_BOTTOM;
    const parts = [];
    const x0 = 4;
    const sysX = single ? x0 + 4 : x0 + 15;
    const clefX = sysX + 8;
    const startX = clefX + (single ? CLEFS.treble.width : Math.max(CLEFS.treble.width, CLEFS.bass.width)) + 16;
    const end = startX + cols * colW + 12;
    const centers = [];

    for (let i = 0; i <= 8; i += 2) {
      parts.push(line('staff-line', sysX, end, GS_TOP + i * H));
      if (!single) parts.push(line('staff-line', sysX, end, GS_B_TOP + i * H));
    }
    if (!single) parts.push(`<path class="brace" d="${bracePath(x0, GS_TOP, GS_B_BOTTOM)}"/>`);
    parts.push(`<line class="barline" x1="${sysX}" x2="${sysX}" y1="${GS_TOP}" y2="${sysBottom}"/>`);
    parts.push(`<line class="barline" x1="${end - 9}" x2="${end - 9}" y1="${GS_TOP}" y2="${sysBottom}"/>`);
    parts.push(`<rect class="barline-thick" x="${end - 5}" y="${GS_TOP}" width="4.5" height="${sysBottom - GS_TOP}"/>`);
    parts.push(glyph('clef', CLEFS.treble, clefX, GS_T_BOTTOM));
    if (!single) parts.push(glyph('clef', CLEFS.bass, clefX, GS_B_BOTTOM));

    // Single staff: room for notes down to C3 and up to B5 with an accidental
    const vbY = single ? GS_T_BOTTOM - 15 * H : 0;
    const vbH = single ? 26 * H + 4 : GS_HEIGHT + (opts.sublabels ? 22 : 0);
    const cursors = [];

    columns.forEach((col, ci) => {
      const center = startX + ci * colW + colW / 2;
      centers[ci] = center;
      if (col && col.cursor) {
        cursors.push(`<rect class="cursor" x="${center - colW / 2 + 3}" y="${vbY + 4}" width="${colW - 6}" height="${vbH - 8}" rx="9"/>`);
      }
      if (col && !col.notes.length && (col.label || col.sublabel)) {
        // labels without notes (e.g. answers typed before the chords are revealed)
        if (col.label) parts.push(`<text class="iv-label ${col.labelState || ''}" x="${center}" y="${GS_LABEL_Y}" text-anchor="middle">${col.label}</text>`);
        return;
      }
      if (!col || !col.notes.length) return;

      parts.push(columnParts(center, col.notes, single));
      if (col.label) {
        parts.push(`<text class="iv-label ${col.labelState || ''}" x="${center}" y="${GS_LABEL_Y}" text-anchor="middle">${col.label}</text>`);
      }
      if (col.sublabel) {
        parts.push(`<text class="iv-label ${col.sublabelState || ''}" x="${center}" y="${GS_LABEL_Y + 22}" text-anchor="middle">${col.sublabel}</text>`);
      }
    });

    svg.innerHTML = cursors.join('') + parts.join('');
    return { scale: sizeSvg(svg, 0, vbY, end, vbH, opts), centers };
  }

  // ---------- Audio: sampled piano, with a synth fallback ----------
  const SAMPLE_BASE = 'https://tonejs.github.io/audio/salamander/';
  let actx = null, master = null, bus = null, ready = Promise.resolve();
  const voices = new Set();    // sources that are playing or scheduled
  const samples = new Map();   // midi -> AudioBuffer

  // Call from a click or key handler. Resolves once the piano samples have loaded.
  function ensureAudio() {
    if (!actx) {
      const AC = window.AudioContext || window.webkitAudioContext;
      if (!AC) return ready;
      actx = new AC();
      master = actx.createGain();
      master.gain.value = 0.9;
      // Keeps big chords from clipping
      const comp = actx.createDynamicsCompressor();
      comp.threshold.value = -10;
      comp.ratio.value = 4;
      master.connect(comp).connect(actx.destination);
      bus = actx.createGain();
      bus.connect(master);
      ready = loadSamples();
    }
    if (actx.state === 'suspended') actx.resume();
    return ready;
  }

  function loadSamples() {
    const names = ['C', 'Ds', 'Fs', 'A'], offs = [0, 3, 6, 9];
    const list = [];
    // Fs1 to C7 covers three ledger lines in every clef, plus accidentals
    for (let oct = 1; oct <= 6; oct++) names.forEach((n, i) => list.push([`${n}${oct}`, 12 * (oct + 1) + offs[i]]));
    list.splice(0, 2);
    list.push(['C7', 96]);
    // The single-file build (build.py) embeds the recordings as base64 in window.EMBEDDED_SAMPLES
    const embedded = window.EMBEDDED_SAMPLES || {};
    const fromBase64 = b64 => Uint8Array.from(atob(b64), c => c.charCodeAt(0)).buffer;
    const load = file => (embedded[file]
      ? Promise.resolve(fromBase64(embedded[file]))
      : fetch(`${SAMPLE_BASE}${file}.mp3`).then(r => { if (!r.ok) throw new Error(r.status); return r.arrayBuffer(); }));
    return Promise.allSettled(list.map(([file, midi]) =>
      load(file)
        .then(buf => new Promise((res, rej) => actx.decodeAudioData(buf, res, rej)))
        .then(audio => samples.set(midi, trim(audio, 3)))));
  }

  // The recordings run ~23 s each but no note rings past ~2.7 s of sample,
  // so keep only the start: about 25 MB of memory instead of about 190 MB.
  function trim(buf, seconds) {
    const len = Math.min(buf.length, Math.floor(seconds * buf.sampleRate));
    const out = actx.createBuffer(buf.numberOfChannels, len, buf.sampleRate);
    for (let c = 0; c < buf.numberOfChannels; c++) out.copyToChannel(buf.getChannelData(c).subarray(0, len), c);
    return out;
  }

  // dur (seconds): damp the note after that long, like lifting the key
  function play(midi, delay = 0, vol = 1, dur = 0) {
    if (!actx || !GAME_SCREENS.has(screen)) return;   // a late timer after leaving stays silent
    const t = actx.currentTime + delay;
    const hold = dur || 1.2, end = dur ? dur + 0.3 : 2.2;
    let nearest = null;
    for (const m of samples.keys()) if (nearest === null || Math.abs(m - midi) < Math.abs(nearest - midi)) nearest = m;

    if (nearest !== null && Math.abs(nearest - midi) <= 3) {
      const src = actx.createBufferSource();
      src.buffer = samples.get(nearest);
      src.playbackRate.value = 2 ** ((midi - nearest) / 12);
      const g = actx.createGain();
      g.gain.setValueAtTime(vol, t);
      g.gain.setValueAtTime(vol, t + hold);
      g.gain.exponentialRampToValueAtTime(0.001, t + end);
      src.connect(g).connect(bus);
      track(src);
      src.start(t);
      src.stop(t + end + 0.1);
      return;
    }
    // Fallback: a few decaying partials
    const f = 440 * 2 ** ((midi - 69) / 12);
    const g = actx.createGain();
    g.gain.setValueAtTime(0, t);
    g.gain.linearRampToValueAtTime(0.45 * vol, t + 0.006);
    g.gain.exponentialRampToValueAtTime(0.001, t + Math.min(end, 1.8));
    g.connect(bus);
    [[1, 0.6], [2, 0.25], [3, 0.1], [4, 0.05]].forEach(([h, amp]) => {
      const o = actx.createOscillator();
      const og = actx.createGain();
      o.frequency.value = f * h;
      og.gain.value = amp;
      o.connect(og).connect(g);
      track(o);
      o.start(t);
      o.stop(t + Math.min(end, 1.8) + 0.1);
    });
  }

  function track(node) {
    voices.add(node);
    node.onended = () => voices.delete(node);
  }

  // Silence everything now, including notes scheduled for later: fade the
  // current bus out quickly (no click), stop its sources, start a fresh bus.
  function stopAll() {
    if (!actx || !bus) return;
    const old = bus, playing = [...voices];
    voices.clear();
    old.gain.setTargetAtTime(0, actx.currentTime, 0.015);
    setTimeout(() => {
      for (const v of playing) { try { v.stop(); } catch (e) { /* already stopped */ } }
      old.disconnect();
    }, 120);
    bus = actx.createGain();
    bus.connect(master);
  }

  // ---------- One-octave piano ----------
  // Builds C–B keys into container; onPress(pitchClass) fires on pointer down
  // (or Enter/Space on a focused key). Returns the key elements by pitch class.
  function makePiano(container, onPress) {
    const keyEls = [];
    const whites = [[0, 'C'], [2, 'D'], [4, 'E'], [5, 'F'], [7, 'G'], [9, 'A'], [11, 'B']];
    const blacks = [[1, 1, 'C sharp / D flat'], [3, 2, 'D sharp / E flat'], [6, 4, 'F sharp / G flat'], [8, 5, 'G sharp / A flat'], [10, 6, 'A sharp / B flat']];
    for (const [pc, name] of whites) {
      const b = document.createElement('button');
      b.type = 'button';
      b.className = 'key white';
      b.textContent = name;
      b.setAttribute('aria-label', name);
      keyEls[pc] = b;
      container.appendChild(b);
    }
    for (const [pc, after, name] of blacks) {
      const b = document.createElement('button');
      b.type = 'button';
      b.className = 'key black';
      b.style.left = `calc(100% / 7 * ${after} - 100% / 7 * 0.3 - 1px)`;
      b.setAttribute('aria-label', name);
      keyEls[pc] = b;
      container.appendChild(b);
    }
    keyEls.forEach((el, pc) => {
      el.addEventListener('pointerdown', e => {
        e.preventDefault();
        el.classList.add('down');
        onPress(pc);
      });
      const up = () => el.classList.remove('down');
      el.addEventListener('pointerup', up);
      el.addEventListener('pointerleave', up);
      el.addEventListener('pointercancel', up);
      el.addEventListener('click', e => { if (e.detail === 0) onPress(pc); });
    });
    return keyEls;
  }

  // Briefly show a key as pressed (for computer-keyboard input)
  function tapKey(el) {
    el.classList.add('down');
    setTimeout(() => el.classList.remove('down'), 120);
  }

  // ---------- Stats (kept in this browser) ----------
  // Per drill: finished rounds (for the progress charts), per-item tallies
  // (notes, intervals, chords) and mix-ups ("answer>guess" counts).
  const STATS_KEY = 'practice.stats.v1';
  const STATS_MODES = ['clefs', 'intervals', 'chords', 'melody', 'progressions', 'intonation', 'standards'];
  let stats = null;

  function readStats() {
    if (!stats) {
      const saved = loadSettings(STATS_KEY, {}, s => typeof s === 'object');
      stats = {};
      for (const m of STATS_MODES) stats[m] = { rounds: [], items: {}, mixups: {}, endless: {}, games: [], ...saved[m] };
    }
    return stats;
  }
  const saveStats = () => saveSettings(STATS_KEY, stats);

  function recordItem(mode, id, ok, ms, meta) {
    const items = readStats()[mode].items;
    const it = items[id] || (items[id] = { n: 0, ok: 0, ms: 0 });
    it.n++;
    if (ok) it.ok++;
    it.ms += Math.round(ms);
    it.meta = meta;
    saveStats();
  }
  function recordMixup(mode, answerId, guessId) {
    const mix = readStats()[mode].mixups;
    const k = `${answerId}>${guessId}`;
    mix[k] = (mix[k] || 0) + 1;
    saveStats();
  }
  function recordRound(mode, round) {
    const rounds = readStats()[mode].rounds;
    rounds.push({ at: Date.now(), ...round });
    if (rounds.length > 500) rounds.shift();
    saveStats();
  }
  // Endless runs, per difficulty ('easy' | 'medium' | 'hard', or 'any' for drills without one)
  function recordEndless(mode, level, run) {
    const all = readStats()[mode].endless;
    const r = all[level] || (all[level] = { best: 0, runs: 0, total: 0 });
    const prev = r.best;
    r.runs++;
    r.total += run;
    r.last = Date.now();
    if (run > r.best) { r.best = run; r.at = Date.now(); }
    saveStats();
    return { best: r.best, prev, isRecord: run > prev };
  }
  // Standard Guesser games: { score, songs, avgMs }. Returns the bests before this game.
  function recordGame(mode, game) {
    const games = readStats()[mode].games;
    const prev = {
      score: Math.max(0, ...games.map(g => g.score)),
      songs: Math.max(0, ...games.map(g => g.songs)),
    };
    games.push({ at: Date.now(), ...game });
    if (games.length > 1000) games.shift();
    saveStats();
    return prev;
  }
  function resetStats() {
    stats = null;
    try { localStorage.removeItem(STATS_KEY); } catch (e) { /* ignore */ }
  }

  // ---------- Results ----------
  let resultsReview = null, resultsAgain = null, resultsSettings = null;

  /**
   * opts: correct, total, times (ms, one per answer), unit, statsKey, summary?, extra?,
   * review (fn(svg, fitPx) or null), onAgain, onSettings,
   * endless?: { level, name? } for an endless run (correct = the run's length)
   * custom?: { title, tiles: [[value, label], [value, label]], summary } shown as given
   */
  function showResults(opts) {
    const { correct, total, times } = opts;
    const avg = times.reduce((a, b) => a + b, 0) / Math.max(1, times.length) / 1000;
    if (opts.custom) {
      const { title, tiles: [[v1, l1], [v2, l2]], summary } = opts.custom;
      $('resultsTitle').textContent = title;
      $('accValue').textContent = v1;
      $('accLabel').textContent = l1;
      $('timeValue').textContent = v2;
      $('timeLabel').textContent = l2;
      $('summary').textContent = summary;
    } else if (opts.endless) {
      const rec = recordEndless(opts.statsKey, opts.endless.level, correct);
      const on = opts.endless.name ? ` on ${opts.endless.name}` : '';
      $('resultsTitle').textContent = 'Run over';
      $('accValue').textContent = String(correct);
      $('accLabel').textContent = 'Run length';
      $('timeValue').textContent = String(rec.best);
      $('timeLabel').textContent = `Longest run${on}`;
      const record = rec.isRecord ? (rec.prev ? `New record${on}, up from ${rec.prev}! ` : `New record${on}! `) : '';
      $('summary').textContent = `${record}${opts.summary || `${correct} ${correct === 1 ? opts.unit : `${opts.unit}s`} right before ${LIVES} misses.`} ${avg.toFixed(2)} s per ${opts.unit}.`;
    } else {
      if (opts.statsKey) recordRound(opts.statsKey, { avgMs: Math.round(avg * 1000), correct, total });
      $('resultsTitle').textContent = 'Round complete';
      $('accValue').textContent = `${Math.round((correct / total) * 100)}%`;
      $('accLabel').textContent = 'Accuracy';
      $('timeValue').textContent = `${avg.toFixed(2)} s`;
      $('timeLabel').textContent = `Avg. time per ${opts.unit}`;
      $('summary').textContent = opts.summary || `${correct} of ${total} ${opts.unit}s correct.`;
    }
    $('extra').textContent = opts.extra || '';
    $('extra').hidden = !opts.extra;
    resultsReview = opts.review || null;
    resultsAgain = opts.onAgain;
    resultsSettings = opts.onSettings;
    $('reviewCard').hidden = !resultsReview;
    show('results');
    drawResultsReview();
  }
  function drawResultsReview() {
    if (resultsReview && screen === 'results') resultsReview($('review'), $('reviewScroll').clientWidth);
  }

  // ---------- Redraw (fonts, resize) ----------
  const redrawers = [drawResultsReview];
  const onRedraw = fn => redrawers.push(fn);
  const redraw = () => redrawers.forEach(fn => fn());

  async function boot() {
    $('againBtn').addEventListener('click', () => resultsAgain && resultsAgain());
    $('settingsBtn').addEventListener('click', () => resultsSettings && resultsSettings());
    document.addEventListener('pointerdown', ensureAudio, { once: true });

    let frame = 0;
    window.addEventListener('resize', () => {
      cancelAnimationFrame(frame);
      frame = requestAnimationFrame(redraw);
    });

    measureGlyphs();
    redraw();
    try {
      const glyphs = [...Object.values(CLEFS), ...Object.values(ACCIDENTALS)].map(g => g.glyph).join('');
      await Promise.race([document.fonts.load('40px "Noto Music"', glyphs), wait(3000)]);
    } catch (e) { /* fall back to system music font */ }
    measureGlyphs();
    redraw();
    // If the web font arrives late, refit the glyphs to it.
    document.fonts.addEventListener?.('loadingdone', () => { measureGlyphs(); redraw(); });
  }

  return {
    LETTERS, SEMIS, CLEFS, ACCIDENTALS, MUSIC_FONT, dia, diaName, diaMidi, accFor, noteName,
    $, randInt, pick, wait, reduceMotion, loadSettings, saveSettings, makeTimers,
    show, current, onShow, setProgress, LIVES, setRun, countPicker,
    renderStaff, renderGrandStaff, roomFor, fitInStage, watchStage,
    ensureAudio, play, stopAll, makePiano, tapKey,
    showResults, onRedraw, boot,
    readStats, recordItem, recordMixup, recordRound, recordEndless, recordGame, resetStats,
    // Drawing pieces for the stats page (staff units: one step = H)
    draw: { H, line, glyph, notePath, ledgers, columnParts, trebleTop: GS_TOP },
  };
})();
