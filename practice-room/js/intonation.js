// Intonation drill: hear a named interval whose top note is slightly out of tune, and say flat or sharp.
(() => {
  'use strict';
  const { diaMidi, $, pick } = App;

  // just: 5-limit ratio; steps: letter names spanned (the tritone as an augmented 4th or a diminished 5th)
  const INTERVALS = [
    { id: 'P1', name: 'Unison',      steps: 0, semis: 0,  just: [1, 1] },
    { id: 'm2', name: 'Minor 2nd',   steps: 1, semis: 1,  just: [16, 15] },
    { id: 'M2', name: 'Major 2nd',   steps: 1, semis: 2,  just: [9, 8] },
    { id: 'm3', name: 'Minor 3rd',   steps: 2, semis: 3,  just: [6, 5] },
    { id: 'M3', name: 'Major 3rd',   steps: 2, semis: 4,  just: [5, 4] },
    { id: 'P4', name: 'Perfect 4th', steps: 3, semis: 5,  just: [4, 3] },
    { id: 'TT', name: 'Tritone',     steps: 3, semis: 6,  just: [45, 32], alt: { steps: 4, just: [64, 45], name: 'Diminished 5th' }, altOf: 'Augmented 4th' },
    { id: 'P5', name: 'Perfect 5th', steps: 4, semis: 7,  just: [3, 2] },
    { id: 'm6', name: 'Minor 6th',   steps: 5, semis: 8,  just: [8, 5] },
    { id: 'M6', name: 'Major 6th',   steps: 5, semis: 9,  just: [5, 3] },
    { id: 'm7', name: 'Minor 7th',   steps: 6, semis: 10, just: [16, 9] },
    { id: 'M7', name: 'Major 7th',   steps: 6, semis: 11, just: [15, 8] },
    { id: 'P8', name: 'Octave',      steps: 7, semis: 12, just: [2, 1] },
  ];
  const BY_ID = Object.fromEntries(INTERVALS.map(iv => [iv.id, iv]));
  const LEVELS = {
    easy:   { name: 'Easy',   cents: [30, 45], desc: '30–45 cents off' },
    medium: { name: 'Medium', cents: [15, 30], desc: '15–30 cents off' },
    hard:   { name: 'Hard',   cents: [6, 15],  desc: '6–15 cents off' },
  };
  const TUNINGS = { equal: 'Equal temperament', just: 'Just intonation' };
  const COUNTS = [10, 20, 30];
  const LOW = 55, HIGH = 81;              // G3 to A5
  const BOTTOM_DUR = 0.9, TOP_AT = 1.0, TOP_DUR = 1.5;
  const KEY = 'intonation.settings';

  const settings = App.loadSettings(KEY, { intervals: ['M3', 'P4', 'P5', 'P8'], level: 'easy', tuning: 'equal', count: 10 },
    s => Array.isArray(s.intervals) && s.level in LEVELS && s.tuning in TUNINGS && COUNTS.includes(s.count));
  settings.intervals = settings.intervals.filter(id => id in BY_ID);
  settings.endless = settings.endless === true;
  const save = () => App.saveSettings(KEY, settings);

  const cents = ([n, d]) => 1200 * Math.log2(n / d);

  // ---------- Setup ----------
  const ivGrid = $('inIvGrid'), levelGrid = $('inLevelGrid'), tuningGroup = $('inTuningGroup'), countGroup = $('inCountGroup');
  const startBtn = $('inStartBtn'), hint = $('inHint'), allBtn = $('inAll');

  const syncCount = App.countPicker(countGroup, COUNTS, settings, () => { save(); syncSetup(); }, 'incount');

  function buildSetup() {
    ivGrid.innerHTML = INTERVALS.map(iv => `
      <button type="button" class="chip iv-chip" id="in-${iv.id}" data-iv="${iv.id}" aria-pressed="false" title="${iv.name}">
        <span class="chip-short">${iv.id}</span>
        <span class="chip-long">${iv.name}</span>
      </button>`).join('');
    levelGrid.innerHTML = Object.entries(LEVELS).map(([id, L]) => `
      <button type="button" class="level-tile" role="radio" id="inlevel-${id}" data-level="${id}" aria-checked="false">
        <span class="level-name">${L.name}</span>
        <span class="level-notes">${L.desc}</span>
      </button>`).join('');
    tuningGroup.innerHTML = Object.entries(TUNINGS).map(([id, name]) =>
      `<button type="button" role="radio" id="tuning-${id}" data-tuning="${id}" aria-checked="false">${name}</button>`).join('');

    ivGrid.addEventListener('click', e => {
      const b = e.target.closest('[data-iv]');
      if (!b) return;
      const on = new Set(settings.intervals);
      on.has(b.dataset.iv) ? on.delete(b.dataset.iv) : on.add(b.dataset.iv);
      settings.intervals = INTERVALS.map(iv => iv.id).filter(id => on.has(id));
      save();
      syncSetup();
    });
    allBtn.addEventListener('click', () => {
      settings.intervals = settings.intervals.length === INTERVALS.length ? [] : INTERVALS.map(iv => iv.id);
      save();
      syncSetup();
    });
    const radio = (group, attr, key, parse = x => x) => group.addEventListener('click', e => {
      const b = e.target.closest(`[data-${attr}]`);
      if (!b) return;
      settings[key] = parse(b.dataset[attr]);
      save();
      syncSetup();
    });
    radio(levelGrid, 'level', 'level');
    radio(tuningGroup, 'tuning', 'tuning');
    radio(countGroup, 'count', 'count', Number);
  }

  function syncSetup() {
    for (const b of ivGrid.querySelectorAll('[data-iv]')) b.setAttribute('aria-pressed', String(settings.intervals.includes(b.dataset.iv)));
    for (const b of levelGrid.querySelectorAll('[data-level]')) b.setAttribute('aria-checked', String(b.dataset.level === settings.level));
    for (const b of tuningGroup.querySelectorAll('[data-tuning]')) b.setAttribute('aria-checked', String(b.dataset.tuning === settings.tuning));
    syncCount();
    allBtn.textContent = settings.intervals.length === INTERVALS.length ? 'Clear all' : 'Select all';
    startBtn.disabled = !settings.intervals.length;
    hint.hidden = !!settings.intervals.length;
  }

  // ---------- Questions ----------
  function makeQuestion(prev) {
    const level = LEVELS[settings.level];
    for (;;) {
      const iv = BY_ID[pick(settings.intervals)];
      const spelled = iv.alt && Math.random() < 0.5 ? { ...iv, ...iv.alt } : { ...iv, name: iv.altOf || iv.name };
      const naturals = [];
      for (let d = 0; d < 70; d++) if (diaMidi(d) >= LOW && diaMidi(d) + iv.semis <= HIGH) naturals.push(d);
      const bd = pick(naturals);
      const bottom = { d: bd, alter: 0, midi: diaMidi(bd) };
      const td = bd + spelled.steps, tm = bottom.midi + iv.semis;
      const top = { d: td, alter: tm - diaMidi(td), midi: tm };
      if (Math.abs(top.alter) > 1) continue;
      if (prev && prev.iv === iv && prev.bottom.midi === bottom.midi) continue;
      const target = settings.tuning === 'just' ? cents(spelled.just) : iv.semis * 100;
      const [lo, hi] = level.cents;
      const off = Math.round(lo + Math.random() * (hi - lo)) * pick([1, -1]);
      return { answerId: off > 0 ? 'sharp' : 'flat', iv, spelled, bottom, top, target, off };
    }
  }

  // Bottom note, then the top note: in tune (off = 0) or detuned by q.off cents
  function playPair(q, off) {
    App.play(q.bottom.midi, 0, 1, BOTTOM_DUR);
    App.play(q.bottom.midi + (q.target + off) / 100, TOP_AT, 1, TOP_DUR);
  }

  const signed = c => `${c > 0 ? '+' : '−'}${Math.abs(c)}¢`;

  function caption(q) {
    const ref = settings.tuning === 'just' ? `just ${q.spelled.just.join(':')}` : 'equal temperament';
    return `${q.spelled.name} (${q.iv.id}) · ${ref}`;
  }

  function columns(q) {
    const answered = !!q.answer;
    const state = answered ? (q.answer.ok ? 'ok' : 'bad') : '';
    return [
      { notes: [{ ...q.bottom, state: '' }] },
      {
        notes: [{ ...q.top, state }],
        label: q.iv.id,
        ...(answered ? { sublabel: `${signed(q.off)} ${q.off > 0 ? 'sharp' : 'flat'}`, sublabelState: state } : {}),
      },
      null,
    ];
  }

  function describeMisses(misses) {
    const sharpAsFlat = misses.filter(([a]) => a === 'sharp').length;
    const flatAsSharp = misses.length - sharpAsFlat;
    if (!misses.length) return '';
    const parts = [];
    if (sharpAsFlat) parts.push(`${sharpAsFlat} sharp ${sharpAsFlat === 1 ? 'note' : 'notes'} as flat`);
    if (flatAsSharp) parts.push(`${flatAsSharp} flat ${flatAsSharp === 1 ? 'note' : 'notes'} as sharp`);
    return `You heard ${parts.join(' and ')}.`;
  }

  function start() {
    EarGame.start({
      total: settings.count,
      unit: 'interval',
      statsKey: 'intonation',
      endless: settings.endless ? { level: settings.level, name: LEVELS[settings.level].name } : null,
      colWidth: 96,
      sublabels: true,
      rightPause: 1600,
      noMixups: true,
      keys: { f: 'flat', s: 'sharp', ArrowLeft: 'flat', ArrowRight: 'sharp' },
      choices: [
        { id: 'flat', label: '♭ Flat', name: 'Flat' },
        { id: 'sharp', label: '♯ Sharp', name: 'Sharp' },
      ],
      statsItem: q => ({ id: q.iv.id, label: q.iv.id }),
      makeQuestion,
      play: q => playPair(q, q.off),
      extra: { label: 'In tune', play: q => playPair(q, 0) },
      caption, columns, describeMisses,
      onAgain: start,
      onSettings: openSetup,
    });
  }

  function openSetup() {
    EarGame.stop();
    syncSetup();
    App.show('intonationSetup');
  }

  buildSetup();
  syncSetup();
  startBtn.addEventListener('click', start);

  window.Intonation = {
    open: openSetup,
    // For the stats page: the same pictures as the Intervals drill
    order: INTERVALS.map(iv => iv.id),
    shape: id => (window.Intervals ? window.Intervals.shape(id) : []),
  };
})();
