// Chords drill: hear a chord, see its first note on a grand staff, and name its quality.
(() => {
  'use strict';
  const { dia, diaMidi, $, pick, randInt } = App;

  // tones: [letter steps above the root, semitones above the root], lowest first.
  // 13 and m11 leave out the 5th, as they usually are voiced.
  const SECTIONS = [
    { id: 'triads', name: 'Triads', chords: [
      { id: 'maj', label: 'maj', name: 'Major triad', short: 'Major', tones: [[0, 0], [2, 4], [4, 7]] },
      { id: 'min', label: 'm', name: 'Minor triad', short: 'Minor', tones: [[0, 0], [2, 3], [4, 7]] },
      { id: 'dim', label: 'dim', name: 'Diminished triad', short: 'Diminished', tones: [[0, 0], [2, 3], [4, 6]] },
      { id: 'aug', label: 'aug', name: 'Augmented triad', short: 'Augmented', tones: [[0, 0], [2, 4], [4, 8]] },
    ] },
    { id: 'sevenths', name: 'Seventh chords', chords: [
      { id: 'dom7', label: '7', name: 'Dominant 7th', short: 'Dominant', tones: [[0, 0], [2, 4], [4, 7], [6, 10]] },
      { id: 'maj7', label: 'maj7', name: 'Major 7th', short: 'Major 7', tones: [[0, 0], [2, 4], [4, 7], [6, 11]] },
      { id: 'min7', label: 'm7', name: 'Minor 7th', short: 'Minor 7', tones: [[0, 0], [2, 3], [4, 7], [6, 10]] },
      { id: 'hdim7', label: 'm7♭5', name: 'Half-diminished 7th', short: 'Half-dim.', tones: [[0, 0], [2, 3], [4, 6], [6, 10]] },
      { id: 'dim7', label: 'dim7', name: 'Diminished 7th', short: 'Dim. 7', tones: [[0, 0], [2, 3], [4, 6], [6, 9]] },
      { id: 'minmaj7', label: 'mMaj7', name: 'Minor-major 7th', short: 'Min-maj 7', tones: [[0, 0], [2, 3], [4, 7], [6, 11]] },
      { id: 'aug7', label: 'aug7', name: 'Augmented 7th', short: 'Aug. 7', tones: [[0, 0], [2, 4], [4, 8], [6, 10]] },
    ] },
    { id: 'extended', name: 'Extensions', chords: [
      { id: 'dom9', label: '9', name: 'Dominant 9th', short: 'Dom. 9', tones: [[0, 0], [2, 4], [4, 7], [6, 10], [8, 14]] },
      { id: 'maj9', label: 'maj9', name: 'Major 9th', short: 'Major 9', tones: [[0, 0], [2, 4], [4, 7], [6, 11], [8, 14]] },
      { id: 'min9', label: 'm9', name: 'Minor 9th', short: 'Minor 9', tones: [[0, 0], [2, 3], [4, 7], [6, 10], [8, 14]] },
      { id: 'dom7b9', label: '7♭9', name: 'Dominant 7th flat 9', short: 'Flat 9', tones: [[0, 0], [2, 4], [4, 7], [6, 10], [8, 13]] },
      { id: 'dom7s9', label: '7♯9', name: 'Dominant 7th sharp 9', short: 'Sharp 9', tones: [[0, 0], [2, 4], [4, 7], [6, 10], [8, 15]] },
      { id: 'min11', label: 'm11', name: 'Minor 11th', short: 'Minor 11', tones: [[0, 0], [2, 3], [6, 10], [8, 14], [10, 17]] },
      { id: 'dom13', label: '13', name: 'Dominant 13th', short: 'Dom. 13', tones: [[0, 0], [2, 4], [6, 10], [8, 14], [12, 21]] },
    ] },
    { id: 'sus6', name: 'Sus & 6 chords', chords: [
      { id: 'sus2', label: 'sus2', name: 'Suspended 2nd', short: 'Sus 2', tones: [[0, 0], [1, 2], [4, 7]] },
      { id: 'sus4', label: 'sus4', name: 'Suspended 4th', short: 'Sus 4', tones: [[0, 0], [3, 5], [4, 7]] },
      { id: 'dom7sus4', label: '7sus4', name: 'Dominant 7th sus 4', short: '7 sus 4', tones: [[0, 0], [3, 5], [4, 7], [6, 10]] },
      { id: 'maj6', label: '6', name: 'Major 6th', short: 'Major 6', tones: [[0, 0], [2, 4], [4, 7], [5, 9]] },
      { id: 'min6', label: 'm6', name: 'Minor 6th', short: 'Minor 6', tones: [[0, 0], [2, 3], [4, 7], [5, 9]] },
      { id: 'six9', label: '6/9', name: 'Six-nine', short: 'Six-nine', tones: [[0, 0], [2, 4], [4, 7], [5, 9], [8, 14]] },
      { id: 'add9', label: 'add9', name: 'Added 9th', short: 'Add 9', tones: [[0, 0], [2, 4], [4, 7], [8, 14]] },
    ] },
  ];
  const CHORDS = SECTIONS.flatMap(s => s.chords);
  const ORDER = CHORDS.map(c => c.id);
  const BY_ID = Object.fromEntries(CHORDS.map(c => [c.id, c]));
  const DIRS = [
    { id: 'asc',  name: 'Ascending' },
    { id: 'desc', name: 'Descending' },
    { id: 'sim',  name: 'Simultaneous' },
  ];
  const ROOTS = [['C', 0], ['D', -1], ['D', 0], ['E', -1], ['E', 0], ['F', 0], ['F', 1], ['G', 0], ['A', -1], ['A', 0], ['B', -1], ['B', 0]];
  const COUNTS = [5, 10, 20, 30];
  const LOW = 41, HIGH = 84;              // F2 to C6
  const GAP = 0.45;                       // seconds between notes of a broken chord
  const OPEN_CHANCE = 0.4;
  const KEY = 'chords.settings';

  const settings = App.loadSettings(KEY, { types: ['maj', 'min', 'dim', 'aug'], dirs: ['sim'], inversions: false, count: 10 },
    s => Array.isArray(s.types) && Array.isArray(s.dirs) && COUNTS.includes(s.count));
  settings.types = settings.types.filter(id => id in BY_ID);
  settings.dirs = settings.dirs.filter(id => DIRS.some(d => d.id === id));
  settings.inversions = settings.inversions === true;
  settings.endless = settings.endless === true;
  const save = () => App.saveSettings(KEY, settings);

  const sectionsEl = $('chordSections'), dirGroup = $('chDirGroup'), countGroup = $('chCountGroup');
  const startBtn = $('chStartBtn'), hint = $('chHint'), allBtn = $('chAll'), invToggle = $('invToggle');

  // ---------- Setup ----------
  function toggle(list, id, order) {
    const on = new Set(list);
    on.has(id) ? on.delete(id) : on.add(id);
    return order.filter(x => on.has(x));
  }

  const syncCount = App.countPicker(countGroup, COUNTS, settings, () => { save(); syncSetup(); }, 'chcount');

  function buildSetup() {
    sectionsEl.innerHTML = SECTIONS.map(s => `
      <div class="chord-section">
        <div class="section-head">
          <h3>${s.name}</h3>
          <button type="button" class="text-btn small" id="sec-${s.id}" data-section="${s.id}">All</button>
        </div>
        <div class="chord-row">
          ${s.chords.map(c => `
            <button type="button" class="chip chord-chip" id="ch-${c.id}" data-chord="${c.id}" aria-pressed="false" title="${c.name}" aria-label="${c.name}">
              <span class="chip-short">${c.label}</span>
              <span class="chip-long">${c.short}</span>
            </button>`).join('')}
        </div>
      </div>`).join('');
    dirGroup.innerHTML = DIRS.map(d =>
      `<button type="button" class="chip dir-chip" id="chdir-${d.id}" data-dir="${d.id}" aria-pressed="false">${d.name}</button>`).join('');

    sectionsEl.addEventListener('click', e => {
      const chip = e.target.closest('[data-chord]');
      const sec = e.target.closest('[data-section]');
      if (chip) {
        settings.types = toggle(settings.types, chip.dataset.chord, ORDER);
      } else if (sec) {
        const ids = SECTIONS.find(s => s.id === sec.dataset.section).chords.map(c => c.id);
        const allOn = ids.every(id => settings.types.includes(id));
        const on = new Set(settings.types);
        ids.forEach(id => (allOn ? on.delete(id) : on.add(id)));
        settings.types = ORDER.filter(id => on.has(id));
      } else return;
      save();
      syncSetup();
    });
    allBtn.addEventListener('click', () => {
      settings.types = settings.types.length === CHORDS.length ? [] : [...ORDER];
      save();
      syncSetup();
    });
    dirGroup.addEventListener('click', e => {
      const b = e.target.closest('[data-dir]');
      if (!b) return;
      settings.dirs = toggle(settings.dirs, b.dataset.dir, DIRS.map(d => d.id));
      save();
      syncSetup();
    });
    invToggle.addEventListener('click', () => {
      settings.inversions = !settings.inversions;
      save();
      syncSetup();
    });
  }

  function syncSetup() {
    for (const b of sectionsEl.querySelectorAll('[data-chord]')) b.setAttribute('aria-pressed', String(settings.types.includes(b.dataset.chord)));
    for (const s of SECTIONS) {
      const allOn = s.chords.every(c => settings.types.includes(c.id));
      $(`sec-${s.id}`).textContent = allOn ? 'None' : 'All';
    }
    for (const b of dirGroup.querySelectorAll('[data-dir]')) b.setAttribute('aria-pressed', String(settings.dirs.includes(b.dataset.dir)));
    syncCount();
    invToggle.setAttribute('aria-checked', String(settings.inversions));
    $('invDesc').textContent = settings.inversions ? 'Inversions & open voicings' : 'Root position only';
    allBtn.textContent = settings.types.length === CHORDS.length ? 'Clear all' : 'Select all';

    const missing = [];
    if (!settings.types.length) missing.push('one chord type');
    if (!settings.dirs.length) missing.push('one direction');
    startBtn.disabled = missing.length > 0;
    hint.hidden = !missing.length;
    hint.textContent = `Pick at least ${missing.join(' and ')} to begin.`;
  }

  // ---------- Voicing ----------
  // Inversions: move the lowest k tones up an octave (extended chords: first inversion at most).
  const maxInversion = n => (n <= 4 ? n - 1 : 1);
  const up = (t, octaves) => ({ ds: t.ds + 7 * octaves, ms: t.ms + 12 * octaves });

  // Relative voicing [{ds, ms}], lowest first. open: spread triad or drop-2 for four-note chords.
  function voice(type, inv, open) {
    let tones = type.tones.map(([ds, ms]) => ({ ds, ms }));
    const k = Math.min(inv, maxInversion(tones.length));
    tones = tones.map((t, i) => (i < k ? up(t, 1) : t)).sort((a, b) => a.ms - b.ms);
    if (open && tones.length === 3) tones[1] = up(tones[1], 1);
    if (open && tones.length === 4) tones[2] = up(tones[2], -1);
    return tones.sort((a, b) => a.ms - b.ms);
  }

  // Absolute notes { d, alter, midi } for a voicing on a root.
  const build = (root, rel) => rel.map(r => {
    const d = root.d + r.ds, midi = root.midi + r.ms;
    return { d, midi, alter: midi - diaMidi(d) };
  });

  // ---------- Questions ----------
  function makeQuestion(prev) {
    for (;;) {
      const type = BY_ID[pick(settings.types)];
      const mode = pick(settings.dirs);
      const n = type.tones.length;
      const inv = settings.inversions ? randInt(0, maxInversion(n)) : 0;
      const open = settings.inversions && n <= 4 && Math.random() < OPEN_CHANCE;
      const rel = voice(type, inv, open);

      const [letter, alter] = pick(ROOTS);
      const roots = [];
      for (let oct = 1; oct <= 5; oct++) {
        const d = dia(letter, oct), midi = diaMidi(d) + alter;
        if (midi + rel[0].ms >= LOW && midi + rel[rel.length - 1].ms <= HIGH) roots.push({ d, midi });
      }
      if (!roots.length) continue;
      const root = pick(roots);
      const notes = build(root, rel);
      if (notes.some(x => Math.abs(x.alter) > 1)) continue;   // skip roots that need double sharps/flats
      if (prev && prev.type === type && prev.root.midi === root.midi) continue;

      // Melodic: the shown note is the first one heard. Together: top or bottom at random.
      const side = mode === 'asc' ? 'bottom' : mode === 'desc' ? 'top' : pick(['bottom', 'top']);
      const shown = side === 'bottom' ? notes[0] : notes[notes.length - 1];
      return { answerId: type.id, type, mode, side, inv, open, root, notes, shown };
    }
  }

  function play(q) {
    if (q.mode === 'sim') {
      const vol = 1.3 / Math.sqrt(q.notes.length);
      for (const x of q.notes) App.play(x.midi, 0, vol);
      return;
    }
    const order = q.mode === 'desc' ? [...q.notes].reverse() : q.notes;
    order.forEach((x, i) => App.play(x.midi, i * GAP));
  }

  function caption(q) {
    const base = q.mode === 'asc' ? 'Ascending' : q.mode === 'desc' ? 'Descending' : 'Simultaneous';
    return q.mode === 'sim' ? `${base} · ${q.side} note shown` : base;
  }

  // Grand-staff columns: shown note, then (once answered) the chord and any wrong guess,
  // the guess built on the same root with the same inversion and spacing.
  function columns(q, guessId) {
    const cols = [{ notes: [{ ...q.shown, state: '' }] }, null, null];
    if (q.answer) {
      cols[1] = { notes: q.notes.map(x => ({ ...x, state: 'ok' })), label: q.type.label, labelState: 'ok' };
    }
    if (guessId) {
      const g = BY_ID[guessId];
      const rel = voice(g, q.inv, q.open && g.tones.length <= 4);
      cols[2] = { notes: build(q.root, rel).map(x => ({ ...x, state: 'bad' })), label: g.label, labelState: 'bad' };
    }
    return cols;
  }

  function start() {
    EarGame.start({
      total: settings.count,
      unit: 'chord',
      statsKey: 'chords',
      endless: settings.endless ? { level: 'any' } : null,
      colWidth: 118,
      choices: settings.types.map(id => ({ id, label: BY_ID[id].label, name: BY_ID[id].name })),
      makeQuestion, play, caption, columns,
      onAgain: start,
      onSettings: openSetup,
    });
  }

  function openSetup() {
    EarGame.stop();
    syncSetup();
    App.show('chordsSetup');
  }

  buildSetup();
  syncSetup();
  startBtn.addEventListener('click', start);

  window.Chords = {
    open: openSetup,
    // For the stats page: chord types in menu order, and each built on middle C
    order: ORDER,
    label: id => BY_ID[id] && BY_ID[id].label,
    shape: id => (BY_ID[id] ? build({ d: 28, midi: 60 }, voice(BY_ID[id], 0, false)) : []),
  };
})();
