// Intervals drill: hear an interval, see its first note on a grand staff, and name it.
(() => {
  'use strict';
  const { diaMidi, $, pick } = App;

  // steps = letter names spanned; the tritone is spelled as an augmented 4th or a diminished 5th.
  const INTERVALS = [
    { id: 'P1', name: 'Unison',      steps: [0], semis: 0 },
    { id: 'm2', name: 'Minor 2nd',   steps: [1], semis: 1 },
    { id: 'M2', name: 'Major 2nd',   steps: [1], semis: 2 },
    { id: 'm3', name: 'Minor 3rd',   steps: [2], semis: 3 },
    { id: 'M3', name: 'Major 3rd',   steps: [2], semis: 4 },
    { id: 'P4', name: 'Perfect 4th', steps: [3], semis: 5 },
    { id: 'TT', name: 'Tritone',     steps: [3, 4], semis: 6 },
    { id: 'P5', name: 'Perfect 5th', steps: [4], semis: 7 },
    { id: 'm6', name: 'Minor 6th',   steps: [5], semis: 8 },
    { id: 'M6', name: 'Major 6th',   steps: [5], semis: 9 },
    { id: 'm7', name: 'Minor 7th',   steps: [6], semis: 10 },
    { id: 'M7', name: 'Major 7th',   steps: [6], semis: 11 },
    { id: 'P8', name: 'Octave',      steps: [7], semis: 12 },
  ];
  const BY_ID = Object.fromEntries(INTERVALS.map(iv => [iv.id, iv]));
  const DIRS = [
    { id: 'asc',  name: 'Ascending' },
    { id: 'desc', name: 'Descending' },
    { id: 'sim',  name: 'Simultaneous' },
  ];
  const COUNTS = [5, 10, 20, 30];
  const LOW = 43, HIGH = 81;              // G2 to A5: every note stays within one ledger line
  const GAP = 0.8;                        // seconds between the two notes of a melodic interval
  const KEY = 'intervals.settings';

  const settings = App.loadSettings(KEY, { intervals: ['m3', 'M3', 'P4', 'P5', 'P8'], dirs: ['asc'], count: 10 },
    s => Array.isArray(s.intervals) && Array.isArray(s.dirs) && COUNTS.includes(s.count));
  settings.intervals = settings.intervals.filter(id => id in BY_ID);
  settings.dirs = settings.dirs.filter(id => DIRS.some(d => d.id === id));
  const save = () => App.saveSettings(KEY, settings);

  const ivGrid = $('ivGrid'), dirGroup = $('dirGroup'), countGroup = $('ivCountGroup');
  const startBtn = $('ivStartBtn'), hint = $('ivHint'), allBtn = $('ivAll');

  // ---------- Setup ----------
  function toggle(list, id, order) {
    const on = new Set(list);
    on.has(id) ? on.delete(id) : on.add(id);
    return order.filter(x => on.has(x));
  }

  function buildSetup() {
    ivGrid.innerHTML = INTERVALS.map(iv => `
      <button type="button" class="chip iv-chip" id="iv-${iv.id}" data-iv="${iv.id}" aria-pressed="false" title="${iv.name}">
        <span class="chip-short">${iv.id}</span>
        <span class="chip-long">${iv.name}</span>
      </button>`).join('');
    dirGroup.innerHTML = DIRS.map(d =>
      `<button type="button" class="chip dir-chip" id="dir-${d.id}" data-dir="${d.id}" aria-pressed="false">${d.name}</button>`).join('');
    countGroup.innerHTML = COUNTS.map(n =>
      `<button type="button" role="radio" id="ivcount-${n}" data-count="${n}" aria-checked="false">${n}</button>`).join('');

    ivGrid.addEventListener('click', e => {
      const b = e.target.closest('[data-iv]');
      if (!b) return;
      settings.intervals = toggle(settings.intervals, b.dataset.iv, INTERVALS.map(iv => iv.id));
      save();
      syncSetup();
    });
    allBtn.addEventListener('click', () => {
      settings.intervals = settings.intervals.length === INTERVALS.length ? [] : INTERVALS.map(iv => iv.id);
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
    countGroup.addEventListener('click', e => {
      const b = e.target.closest('[data-count]');
      if (!b) return;
      settings.count = Number(b.dataset.count);
      save();
      syncSetup();
    });
  }

  function syncSetup() {
    for (const b of ivGrid.querySelectorAll('[data-iv]')) b.setAttribute('aria-pressed', String(settings.intervals.includes(b.dataset.iv)));
    for (const b of dirGroup.querySelectorAll('[data-dir]')) b.setAttribute('aria-pressed', String(settings.dirs.includes(b.dataset.dir)));
    for (const b of countGroup.querySelectorAll('[data-count]')) b.setAttribute('aria-checked', String(Number(b.dataset.count) === settings.count));
    allBtn.textContent = settings.intervals.length === INTERVALS.length ? 'Clear all' : 'Select all';

    const missing = [];
    if (!settings.intervals.length) missing.push('one interval');
    if (!settings.dirs.length) missing.push('one direction');
    startBtn.disabled = missing.length > 0;
    hint.hidden = !missing.length;
    hint.textContent = `Pick at least ${missing.join(' and ')} to begin.`;
  }

  // ---------- Questions ----------
  // A note is { d, alter, midi }: diatonic number, sharps(+)/flats(-), pitch.
  function spell(from, iv, dir, steps) {
    const d = from.d + dir * steps;
    const midi = from.midi + dir * iv.semis;
    return { d, alter: midi - diaMidi(d), midi };
  }

  // The shown note is always natural and has room for an octave in its direction,
  // so both the answer and any wrong guess fit on the staff.
  function makeQuestion(prev) {
    for (;;) {
      const iv = BY_ID[pick(settings.intervals)];
      const mode = pick(settings.dirs);
      // dir: which way the hidden note lies from the shown one
      const dir = mode === 'asc' ? 1 : mode === 'desc' ? -1 : pick([1, -1]);
      const lo = dir === 1 ? LOW : LOW + 12, hi = dir === 1 ? HIGH - 12 : HIGH;
      const naturals = [];
      for (let d = 0; d < 70; d++) if (diaMidi(d) >= lo && diaMidi(d) <= hi) naturals.push(d);
      const sd = pick(naturals);
      const shown = { d: sd, alter: 0, midi: diaMidi(sd) };
      const other = spell(shown, iv, dir, pick(iv.steps));
      if (Math.abs(other.alter) > 1) continue;
      if (prev && prev.iv === iv && prev.shown.midi === shown.midi && prev.mode === mode) continue;
      return { answerId: iv.id, iv, mode, dir, shown, other };
    }
  }

  // The shown note is always the one heard first in melodic modes.
  function play(q) {
    App.play(q.shown.midi);
    App.play(q.other.midi, q.mode === 'sim' ? 0 : GAP);
  }

  function caption(q) {
    if (q.mode === 'asc') return 'Ascending';
    if (q.mode === 'desc') return 'Descending';
    return q.dir === 1 ? 'Simultaneous · bottom note shown' : 'Simultaneous · top note shown';
  }

  // Grand-staff columns: shown note, then (once answered) the answer and any wrong guess.
  function columns(q, guessId) {
    const answered = !!q.answer;
    const shown = { ...q.shown, state: '' };
    const other = answered && { ...q.other, state: 'ok' };
    let guess = null;
    if (guessId) {
      const iv = BY_ID[guessId];
      guess = { ...spell(q.shown, iv, q.dir, iv.steps[0]), state: 'bad' };
    }
    const okLabel = answered ? { label: q.iv.id, labelState: 'ok' } : {};
    if (q.mode === 'sim') {
      return [
        { notes: other ? [shown, other] : [shown], ...okLabel },
        guess && { notes: [shown, guess], label: guessId, labelState: 'bad' },
        null,
      ];
    }
    return [
      { notes: [shown] },
      other && { notes: [other], ...okLabel },
      guess && { notes: [guess], label: guessId, labelState: 'bad' },
    ];
  }

  function start() {
    EarGame.start({
      total: settings.count,
      unit: 'interval',
      statsKey: 'intervals',
      colWidth: 84,
      choices: settings.intervals.map(id => ({ id, label: id, name: BY_ID[id].name })),
      makeQuestion, play, caption, columns,
      onAgain: start,
      onSettings: openSetup,
    });
  }

  function openSetup() {
    EarGame.stop();
    syncSetup();
    App.show('intervalsSetup');
  }

  buildSetup();
  syncSetup();
  startBtn.addEventListener('click', start);

  window.Intervals = {
    open: openSetup,
    // For the stats page: every interval in size order, and its notes above middle C
    order: INTERVALS.map(iv => iv.id),
    shape: id => {
      const iv = BY_ID[id], c4 = { d: 28, alter: 0, midi: 60 };
      return iv ? [c4, spell(c4, iv, 1, iv.steps[0])] : [];
    },
  };
})();
