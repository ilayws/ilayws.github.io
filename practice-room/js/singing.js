// Interval Singing: hear the first note of an interval, sing the second into the microphone,
// and see how many cents off you were.
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
    { id: 'asc',  name: 'Ascending',  short: '↑ Up' },
    { id: 'desc', name: 'Descending', short: '↓ Down' },
  ];
  // tol: how many cents off still counts as right
  const LEVELS = {
    easy:   { name: 'Easy',   tol: 50, desc: 'The right note, ±50¢' },
    medium: { name: 'Medium', tol: 25, desc: 'Within ±25¢' },
    hard:   { name: 'Hard',   tol: 10, desc: 'Within ±10¢' },
  };
  // Both notes of every interval stay in this range
  const VOICES = {
    low:  { name: 'Low',    range: [43, 62] },   // G2–D4
    mid:  { name: 'Middle', range: [48, 67] },   // C3–G4
    high: { name: 'High',   range: [55, 74] },   // G3–D5
  };
  const COUNTS = [5, 10, 20, 30];
  const NOTE_DUR = 1.0;                  // seconds the first note sounds
  const LISTEN_AFTER = 1400;             // ms: start listening once the note has died away
  const GAP = 0.8;                       // seconds between the notes of "Hear both"
  const RIGHT_PAUSE = 2200;
  const KEY = 'singing.settings';

  const settings = App.loadSettings(KEY, { intervals: ['m3', 'M3', 'P4', 'P5', 'P8'], dirs: ['asc'], level: 'easy', voice: 'mid', count: 10 },
    s => Array.isArray(s.intervals) && Array.isArray(s.dirs) && s.level in LEVELS && s.voice in VOICES && COUNTS.includes(s.count));
  settings.intervals = settings.intervals.filter(id => id in BY_ID);
  settings.dirs = settings.dirs.filter(id => DIRS.some(d => d.id === id));
  settings.endless = settings.endless === true;
  const save = () => App.saveSettings(KEY, settings);

  // ---------- Setup ----------
  const ivGrid = $('siIvGrid'), dirGroup = $('siDirGroup'), voiceGroup = $('siVoiceGroup'), levelGrid = $('siLevelGrid');
  const startBtn = $('siStartBtn'), hint = $('siHint'), allBtn = $('siAll');

  const syncCount = App.countPicker($('siCountGroup'), COUNTS, settings, () => { save(); syncSetup(); }, 'sicount');

  function toggle(list, id, order) {
    const on = new Set(list);
    on.has(id) ? on.delete(id) : on.add(id);
    return order.filter(x => on.has(x));
  }

  function buildSetup() {
    ivGrid.innerHTML = INTERVALS.map(iv => `
      <button type="button" class="chip iv-chip" id="si-${iv.id}" data-iv="${iv.id}" aria-pressed="false" title="${iv.name}">
        <span class="chip-short">${iv.id}</span>
        <span class="chip-long">${iv.name}</span>
      </button>`).join('');
    dirGroup.innerHTML = DIRS.map(d =>
      `<button type="button" class="chip dir-chip" id="sidir-${d.id}" data-dir="${d.id}" aria-pressed="false" aria-label="${d.name}"><span class="long">${d.name}</span><span class="short" aria-hidden="true">${d.short}</span></button>`).join('');
    voiceGroup.innerHTML = Object.entries(VOICES).map(([id, v]) => {
      const range = `${App.diaName(lowestD(v.range[0]))}–${App.diaName(lowestD(v.range[1]))}`;
      return `<button type="button" role="radio" id="voice-${id}" data-voice="${id}" aria-checked="false" title="${range}" aria-label="${v.name} (${range})">${v.name}</button>`;
    }).join('');
    levelGrid.innerHTML = Object.entries(LEVELS).map(([id, L]) => `
      <button type="button" class="level-tile" role="radio" id="silevel-${id}" data-level="${id}" aria-checked="false">
        <span class="level-name">${L.name}</span>
        <span class="level-notes">${L.desc}</span>
      </button>`).join('');

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
    const radio = (group, attr, key) => group.addEventListener('click', e => {
      const b = e.target.closest(`[data-${attr}]`);
      if (!b) return;
      settings[key] = b.dataset[attr];
      save();
      syncSetup();
    });
    radio(voiceGroup, 'voice', 'voice');
    radio(levelGrid, 'level', 'level');
  }

  // The natural note at or below a pitch (for naming the voice ranges)
  function lowestD(midi) {
    let d = Math.floor((midi / 12) * 7);
    while (diaMidi(d + 1) <= midi) d++;
    while (diaMidi(d) > midi) d--;
    return d;
  }

  function syncSetup(problem) {
    for (const b of ivGrid.querySelectorAll('[data-iv]')) b.setAttribute('aria-pressed', String(settings.intervals.includes(b.dataset.iv)));
    for (const b of dirGroup.querySelectorAll('[data-dir]')) b.setAttribute('aria-pressed', String(settings.dirs.includes(b.dataset.dir)));
    for (const b of voiceGroup.querySelectorAll('[data-voice]')) b.setAttribute('aria-checked', String(b.dataset.voice === settings.voice));
    for (const b of levelGrid.querySelectorAll('[data-level]')) b.setAttribute('aria-checked', String(b.dataset.level === settings.level));
    syncCount();
    allBtn.textContent = settings.intervals.length === INTERVALS.length ? 'Clear all' : 'Select all';

    const missing = [];
    if (!settings.intervals.length) missing.push('one interval');
    if (!settings.dirs.length) missing.push('one direction');
    startBtn.disabled = missing.length > 0;
    hint.textContent = missing.length ? `Pick at least ${missing.join(' and ')} to begin.` : problem || '';
    hint.hidden = !hint.textContent;
  }

  // ---------- Questions ----------
  // A note is { d, alter, midi }: diatonic number, sharps(+)/flats(-), pitch.
  function spell(from, iv, dir, steps) {
    const d = from.d + dir * steps;
    const midi = from.midi + dir * iv.semis;
    return { d, alter: midi - diaMidi(d), midi };
  }

  // The first note is always natural; both notes fit the voice range.
  function makeQuestion(prev) {
    const [lo, hi] = VOICES[settings.voice].range;
    for (;;) {
      const iv = BY_ID[pick(settings.intervals)];
      const mode = pick(settings.dirs);
      const dir = mode === 'asc' ? 1 : -1;
      const naturals = [];
      for (let d = 0; d < 70; d++) {
        const m = diaMidi(d), t = m + dir * iv.semis;
        if (m >= lo && m <= hi && t >= lo && t <= hi) naturals.push(d);
      }
      const fd = pick(naturals);
      const first = { d: fd, alter: 0, midi: diaMidi(fd) };
      const target = spell(first, iv, dir, pick(iv.steps));
      if (Math.abs(target.alter) > 1) continue;
      if (prev && prev.iv === iv && prev.first.midi === first.midi && prev.mode === mode) continue;
      return { iv, mode, dir, first, target, start: 0, result: null, plays: 0 };
    }
  }

  // How a sung pitch class is written: [letter index, alteration]
  const SPELL = [[0, 0], [0, 1], [1, 0], [2, -1], [2, 0], [3, 0], [3, 1], [4, 0], [5, -1], [5, 0], [6, -1], [6, 0]];
  function spellMidi(midi) {
    const [letter, alter] = SPELL[((midi % 12) + 12) % 12];
    return { d: (Math.floor(midi / 12) - 1) * 7 + letter, alter, midi };
  }

  const signed = c => `${c > 0 ? '+' : c < 0 ? '−' : '±'}${Math.abs(Math.round(c))}¢`;
  const leaning = c => (Math.abs(c) <= App.IN_TUNE ? 'in tune' : c > 0 ? 'sharp' : 'flat');

  // ---------- Microphone ----------
  const mic = { ctx: null, stream: null, src: null, analyser: null, buf: null };

  async function openMic() {
    if (mic.stream) return;
    if (!navigator.mediaDevices || !navigator.mediaDevices.getUserMedia) throw Object.assign(new Error('no mic'), { name: 'NotSupportedError' });
    // Raw input: echo cancelling and noise suppression bend the pitch
    const stream = await navigator.mediaDevices.getUserMedia({ audio: { echoCancellation: false, noiseSuppression: false, autoGainControl: false } });
    const AC = window.AudioContext || window.webkitAudioContext;
    if (!mic.ctx) mic.ctx = new AC();
    if (mic.ctx.state === 'suspended') mic.ctx.resume();
    mic.src = mic.ctx.createMediaStreamSource(stream);
    mic.analyser = mic.ctx.createAnalyser();
    mic.analyser.fftSize = mic.ctx.sampleRate > 50000 ? 4096 : 2048;
    mic.buf = new Float32Array(mic.analyser.fftSize);
    mic.src.connect(mic.analyser);
    mic.stream = stream;
  }

  function closeMic() {
    stopListening();
    if (!mic.stream) return;
    mic.stream.getTracks().forEach(t => t.stop());
    mic.src.disconnect();
    mic.stream = null;
  }

  // ---------- Pitch detection (YIN) ----------
  // Returns the frame's loudness and, when it's a clear pitch, its fractional MIDI number.
  const MIN_HZ = 65, MAX_HZ = 1100, YIN_THRESHOLD = 0.15, GATE = 0.006;
  let yinD = null;
  function detect(buf, sr) {
    let rms = 0;
    for (let i = 0; i < buf.length; i++) rms += buf[i] * buf[i];
    rms = Math.sqrt(rms / buf.length);
    if (rms < GATE) return { rms, midi: null };
    const W = buf.length >> 1;
    const maxTau = Math.min(W - 2, Math.floor(sr / MIN_HZ)), minTau = Math.max(2, Math.floor(sr / MAX_HZ));
    if (!yinD || yinD.length < maxTau + 2) yinD = new Float32Array(maxTau + 2);
    const d = yinD;
    d[0] = 1;
    let total = 0;
    for (let tau = 1; tau <= maxTau + 1; tau++) {
      let sum = 0;
      for (let i = 0; i < W; i++) { const x = buf[i] - buf[i + tau]; sum += x * x; }
      total += sum;
      d[tau] = total ? (sum * tau) / total : 1;
    }
    let tau = -1;
    for (let t = minTau; t <= maxTau; t++) {
      if (d[t] < YIN_THRESHOLD) {
        while (t + 1 <= maxTau && d[t + 1] < d[t]) t++;
        tau = t;
        break;
      }
    }
    if (tau < 0) return { rms, midi: null };
    const a = d[tau - 1], b = d[tau], c = d[tau + 1], curve = a + c - 2 * b;
    const exact = tau + (curve ? (a - c) / (2 * curve) : 0);
    return { rms, midi: 69 + 12 * Math.log2(sr / exact / 440) };
  }

  const median = arr => {
    const s = [...arr].sort((x, y) => x - y), m = s.length >> 1;
    return s.length % 2 ? s[m] : (s[m - 1] + s[m]) / 2;
  };
  // The mean of the middle 80%: steady through vibrato, deaf to the odd stray frame
  const trimmedMean = arr => {
    const s = [...arr].sort((x, y) => x - y), cut = Math.floor(s.length * 0.1), mid = s.slice(cut, s.length - cut);
    return mid.reduce((a, b) => a + b, 0) / mid.length;
  };

  // ---------- Listening ----------
  // A sung note counts once it has held steady for HOLD_MS; its pitch is the average
  // after the first ATTACK_MS (where singers scoop into the note).
  const FRAME_MS = 30, HOLD_MS = 700, ATTACK_MS = 150, GAP_MS = 160, JUMP = 0.7;
  let loop = 0, hold = null;

  function listen() {
    stopListening();
    hold = { run: [], lastVoiced: 0 };
    loop = setInterval(tick, FRAME_MS);
  }
  function stopListening() {
    clearInterval(loop);
    loop = 0;
    levelEl.style.width = '0%';
    holdEl.style.width = '0%';
  }

  function tick() {
    if (!mic.analyser || !game || game.state !== 'listen') return stopListening();
    mic.analyser.getFloatTimeDomainData(mic.buf);
    const { rms, midi } = detect(mic.buf, mic.ctx.sampleRate);
    const now = performance.now(), run = hold.run;
    levelEl.style.width = `${Math.min(100, Math.sqrt(rms / 0.25) * 100)}%`;
    if (midi == null) {
      if (run.length && now - hold.lastVoiced > GAP_MS) run.length = 0;
    } else {
      let m = midi;
      if (run.length) {
        const ref = median(run.slice(-5).map(f => f.m));
        m += 12 * Math.round((ref - m) / 12);   // the detector slipping an octave mid-note
        if (Math.abs(m - ref) > JUMP) { run.length = 0; m = midi; }
      }
      run.push({ t: now, m });
      hold.lastVoiced = now;
    }
    const held = run.length ? now - run[0].t : 0;
    holdEl.style.width = `${Math.min(1, held / HOLD_MS) * 100}%`;
    panel.dataset.state = held ? 'hold' : 'listen';
    statusEl.textContent = held ? 'Hold it…' : singPrompt(game.q);
    if (held >= HOLD_MS) grade(trimmedMean(run.filter(f => f.t - run[0].t >= ATTACK_MS).map(f => f.m)));
  }

  // ---------- Game ----------
  // state: 'loading' | 'play' (first note sounding) | 'listen' | 'right' | 'missed' (waiting for Continue)
  const stage = $('siStage'), card = $('siCard'), grandSvg = $('siGrand'), grandWrap = $('siGrandWrap');
  const caption = $('siCaption'), replayBtn = $('siReplay'), bothBtn = $('siBoth'), nextBtn = $('siNext'), stageRow = $('siStageRow');
  const panel = $('siPanel'), statusEl = $('siStatus'), centsEl = $('siCents'), skipBtn = $('siSkip');
  const zoneEl = $('siZone'), markEl = $('siMark'), levelEl = $('siLevel'), holdEl = $('siHold');
  const timers = App.makeTimers();
  let game = null;

  const singPrompt = q => (q.iv.semis === 0 ? 'Sing the same note' : `Sing the ${q.dir > 0 ? 'higher' : 'lower'} note`);

  function showProgress() {
    const fill = $('siProgressFill'), text = $('siProgressText');
    if (game.endless) App.setRun(fill, text, game.lives, game.results.filter(Boolean).length);
    else App.setProgress(fill, text, game.i, game.count);
  }

  function columns(q) {
    const r = q.result;
    if (!r) return [{ notes: [{ ...q.first, state: '' }] }, null, null];
    const state = r.ok ? 'ok' : 'bad';
    const target = {
      notes: [{ ...q.target, state: r.ok ? 'ok' : r.skipped ? '' : 'bad' }],
      label: q.iv.id,
      labelState: r.skipped ? '' : state,
      ...(r.near ? { sublabel: signed(r.off), sublabelState: state } : {}),
    };
    const sung = r.sung && !r.near && { notes: [{ ...r.sung, state: 'bad' }], label: 'You', labelState: 'bad' };
    return [{ notes: [{ ...q.first, state: '' }] }, target, sung || null];
  }

  function draw() {
    if (!game || !game.q) return;
    App.fitInStage(stage, card, grandWrap, maxHeight =>
      App.renderGrandStaff(grandSvg, columns(game.q), { cols: 3, colWidth: 84, fit: grandWrap.clientWidth - 16, min: 0.5, max: 1.7, sublabels: true, maxHeight }));
  }

  function setNext(visible) {
    nextBtn.hidden = !visible;
    nextBtn.textContent = isLast() ? 'See results' : 'Continue';
    stageRow.classList.toggle('has-continue', visible);
  }

  async function start() {
    const loading = App.ensureAudio();
    startBtn.disabled = true;
    try {
      await openMic();
    } catch (e) {
      startBtn.disabled = false;
      const blocked = e && (e.name === 'NotAllowedError' || e.name === 'SecurityError');
      const none = e && (e.name === 'NotFoundError' || e.name === 'OverconstrainedError');
      syncSetup(blocked ? 'The microphone is blocked. Allow this site to use it (in the address bar), then press Start again.'
        : none ? 'No microphone was found. Connect one and press Start again.'
        : 'This browser can’t use the microphone here. Try the online Practice Room in Chrome, Firefox or Safari.');
      if (App.current() === 'results') App.show('singSetup');
      return;
    }
    startBtn.disabled = false;
    timers.clear();
    game = {
      endless: settings.endless,
      level: settings.level,
      count: settings.count,
      lives: App.LIVES,
      i: 0,
      q: null,
      state: 'loading',
      results: [],
      times: [],
      sungOffs: [],
    };
    setNext(false);
    bothBtn.hidden = true;
    App.show('singGame');
    showProgress();
    caption.textContent = 'Loading piano…';
    setPanel('wait', '');
    const g = game;
    await Promise.race([loading, App.wait(5000)]);
    if (game !== g) return;   // left the screen while loading
    nextQuestion();
  }

  function setPanel(state, status, cents = '') {
    panel.dataset.state = state;
    statusEl.textContent = status;
    centsEl.textContent = cents;
  }

  function nextQuestion() {
    game.q = makeQuestion(game.q);
    game.state = 'play';
    setNext(false);
    bothBtn.hidden = true;
    skipBtn.hidden = false;
    caption.textContent = `${DIRS.find(d => d.id === game.q.mode).name} ${game.q.iv.name.toLowerCase()}`;
    if (game.q.iv.semis === 0) caption.textContent = 'Unison';
    showProgress();
    draw();
    setPanel('wait', 'Listen…');
    timers.later(() => game && playFirst(), 250);
  }

  // Play the first note, then listen once it has died away
  function playFirst() {
    const q = game.q, n = ++q.plays;
    stopListening();
    App.stopAll();
    App.play(q.first.midi, 0, 1, NOTE_DUR);
    if (game.state !== 'play' && game.state !== 'listen') return;
    game.state = 'play';
    setPanel('wait', 'Listen…');
    timers.later(() => {
      if (!game || game.q !== q || q.plays !== n || game.state !== 'play') return;
      game.state = 'listen';
      if (!q.start) q.start = performance.now();
      setPanel('listen', singPrompt(q));
      listen();
    }, LISTEN_AFTER);
  }

  function playBoth() {
    const q = game && game.q;
    if (!q || !q.result) return;
    App.stopAll();
    App.play(q.first.midi);
    App.play(q.target.midi, GAP);
  }

  function grade(sung) {
    stopListening();
    const q = game.q, tol = LEVELS[game.level].tol;
    const raw = (sung - q.target.midi) * 100;
    const off = raw - 1200 * Math.round(raw / 1200);   // any octave counts
    const near = Math.abs(off) <= 50;                  // the right note, however out of tune
    const ok = Math.abs(off) <= tol;
    const semis = Math.round(off / 100);
    q.result = { ok, off, near, sung: near ? null : spellMidi(q.target.midi + semis) };
    record(ok, near ? off : null);

    if (near) {
      game.sungOffs.push(off);
      const pos = Math.max(0, Math.min(100, 50 + off));
      markEl.style.left = `${pos}%`;
      zoneEl.style.left = `${50 - tol}%`;
      zoneEl.style.width = `${2 * tol}%`;
      const how = leaning(off);
      const verdict = ok ? (how === 'in tune' ? 'In tune' : `Right, a little ${how}`) : `Too ${how} for ${LEVELS[game.level].name.toLowerCase()} (±${tol}¢)`;
      setPanel(ok ? 'ok' : 'bad', verdict, signed(off));
    } else {
      const name = App.noteName(q.result.sung.d, q.result.sung.alter);
      const k = Math.abs(semis);
      setPanel('wrong', `You sang ${name}: ${k} ${k === 1 ? 'semitone' : 'semitones'} too ${semis > 0 ? 'high' : 'low'}`);
    }
    afterAnswer(ok);
  }

  function skip() {
    if (!game || (game.state !== 'play' && game.state !== 'listen')) return;
    stopListening();
    timers.clear();
    const q = game.q;
    if (!q.start) q.start = performance.now();
    q.result = { ok: false, skipped: true };
    record(false, null);
    setPanel('wrong', `Skipped. It was ${App.noteName(q.target.d, q.target.alter)}`);
    afterAnswer(false);
  }

  function record(ok, off) {
    const q = game.q, ms = performance.now() - q.start;
    game.results.push(ok);
    game.times.push(ms);
    App.recordItem('singing', q.iv.id, ok, ms, { label: q.iv.id });
    if (off != null) App.recordTuning('singing', String(q.target.midi % 12), off);
  }

  function afterAnswer(ok) {
    const q = game.q;
    skipBtn.hidden = true;
    bothBtn.hidden = false;
    draw();
    if (ok) {
      game.state = 'right';
      showProgress();
      timers.later(advance, RIGHT_PAUSE);
      return;
    }
    game.lives--;
    game.state = 'missed';
    showProgress();
    timers.later(() => { if (game && game.q === q) { App.stopAll(); App.play(q.target.midi); } }, 400);   // the right note
    setNext(true);
    nextBtn.focus({ preventScroll: true });
  }

  const isLast = () => !!game && (game.endless ? game.lives <= 0 : game.i + 1 >= game.count);

  function advance() {
    if (!game) return;
    if (isLast()) return finish();
    game.i++;
    nextQuestion();
  }

  function continueAfterMiss() {
    if (!game || game.state !== 'missed') return;
    timers.clear();
    advance();
  }

  function summary(offs) {
    if (!offs.length) return '';
    const avg = offs.reduce((a, b) => a + b, 0) / offs.length;
    const spread = offs.reduce((a, b) => a + Math.abs(b), 0) / offs.length;
    const how = leaning(avg);
    return `On the notes you hit you were ${how === 'in tune' ? 'in tune' : `${Math.abs(Math.round(avg))}¢ ${how}`} on average (${signed(avg)}), ${Math.round(spread)}¢ off either way.`;
  }

  function finish() {
    const done = game;
    timers.clear();
    stopListening();
    game = null;
    App.showResults({
      correct: done.results.filter(Boolean).length,
      total: done.results.length,
      times: done.times,
      unit: 'interval',
      statsKey: 'singing',
      endless: done.endless ? { level: done.level, name: LEVELS[done.level].name } : null,
      extra: summary(done.sungOffs),
      review: null,
      onAgain: start,
      onSettings: openSetup,
    });
  }

  function openSetup() {
    timers.clear();
    game = null;
    closeMic();
    syncSetup();
    App.show('singSetup');
  }

  function replay() {
    App.ensureAudio();
    if (game && game.q && game.state !== 'loading') playFirst();
  }

  // ---------- Wiring ----------
  buildSetup();
  syncSetup();
  startBtn.addEventListener('click', start);
  $('siQuitBtn').addEventListener('click', openSetup);
  replayBtn.addEventListener('click', replay);
  bothBtn.addEventListener('click', playBoth);
  nextBtn.addEventListener('click', continueAfterMiss);
  skipBtn.addEventListener('click', skip);

  // Stop listening when the game is left, and let go of the microphone unless a new round may follow
  App.onShow(name => {
    if (name !== 'singGame') stopListening();
    if (name !== 'singGame' && name !== 'results') closeMic();
  });

  document.addEventListener('keydown', e => {
    if (e.repeat || e.metaKey || e.ctrlKey || e.altKey) return;
    if (App.current() !== 'singGame') return;
    const k = e.key.toLowerCase();
    if (k === 'escape') return openSetup();
    if (k === 'enter') { e.preventDefault(); return continueAfterMiss(); }
    if (k === ' ') { e.preventDefault(); return replay(); }
    if (k === 's') { e.preventDefault(); return skip(); }
  });

  App.watchStage(stage, () => { if (App.current() === 'singGame') draw(); });
  App.onRedraw(() => { if (App.current() === 'singGame') draw(); });

  window.Singing = {
    open: openSetup,
    // For the stats page: the same pictures as the Intervals drill
    order: INTERVALS.map(iv => iv.id),
    shape: id => (window.Intervals ? window.Intervals.shape(id) : []),
  };
})();
