// Melody drill: hear a short melody, then play it back on a one-octave piano.
(() => {
  'use strict';
  const { LETTERS, SEMIS, dia, diaMidi, accFor, $, pick, reduceMotion } = App;

  // ---------- Keys ----------
  const MAJOR = [0, 2, 4, 5, 7, 9, 11], MINOR = [0, 2, 3, 5, 7, 8, 10];
  // [tonic letter, alter, mode, signature (+ sharps / - flats)]
  const KEYS = [
    ['C', 0, 'major', 0], ['G', 0, 'major', 1], ['F', 0, 'major', -1], ['D', 0, 'major', 2], ['B', -1, 'major', -2],
    ['A', 0, 'major', 3], ['E', -1, 'major', -3], ['E', 0, 'major', 4], ['A', -1, 'major', -4],
    ['A', 0, 'minor', 0], ['E', 0, 'minor', 1], ['D', 0, 'minor', -1], ['B', 0, 'minor', 2], ['G', 0, 'minor', -2],
    ['F', 1, 'minor', 3], ['C', 0, 'minor', -3], ['C', 1, 'minor', 4], ['F', 0, 'minor', -4],
  ].map(([letter, alter, mode, sig]) => ({
    letter, alter, mode, sig,
    scale: mode === 'major' ? MAJOR : MINOR,
    name: `${letter}${alter === 1 ? '♯' : alter === -1 ? '♭' : ''} ${mode}`,
  }));

  // ---------- Difficulty ----------
  // leaps: weight for each melodic interval size, in scale steps (0 = repeat, 7 = octave).
  // low/high: allowed range in scale steps around the tonic; span: widest range in one melody.
  const LEVELS = {
    easy: {
      name: 'Easy', notes: 5, desc: 'Major keys · steps and small leaps',
      keys: k => k.mode === 'major' && Math.abs(k.sig) <= 2,
      leaps: { 0: 0.3, 1: 5, 2: 2, 3: 0.4 }, low: -3, high: 5, span: 5, arch: 2.5,
      beat: 0.7, tritone: false, chromatic: 0,
    },
    medium: {
      name: 'Medium', notes: 8, desc: 'Major and minor · leaps up to a 5th',
      keys: k => Math.abs(k.sig) <= (k.mode === 'major' ? 3 : 2),
      leaps: { 0: 0.4, 1: 5, 2: 2.5, 3: 1, 4: 0.8 }, low: -4, high: 7, span: 7, arch: 3.5,
      beat: 0.6, tritone: false, chromatic: 0,
    },
    hard: {
      name: 'Hard', notes: 12, desc: 'Any key · wide leaps and chromatic notes',
      keys: () => true,
      leaps: { 0: 0.3, 1: 5, 2: 2.5, 3: 1.3, 4: 1, 5: 0.5, 7: 0.4 }, low: -5, high: 9, span: 9, arch: 4.5,
      beat: 0.5, tritone: true, chromatic: 0.6,
    },
  };
  const COUNTS = [3, 5, 10];
  const LOW = 55, HIGH = 83;             // G3 to B5, inside the keyboard's octave range
  const OCT_MIN = 3, OCT_MAX = 5;
  const KEY = 'melody.settings';

  const settings = App.loadSettings(KEY, { level: 'easy', count: 5, keyChord: true },
    s => s.level in LEVELS && COUNTS.includes(s.count));
  settings.keyChord = settings.keyChord !== false;
  const save = () => App.saveSettings(KEY, settings);

  // ---------- Melody generator ----------
  // Works in scale steps (0 = tonic). Each melody follows an implied chord progression:
  // strong beats favour chord tones, other notes move by step (passing / neighbour notes),
  // leaps are followed by a step back, the line follows an overall contour, and it ends
  // with a stepwise cadence onto the tonic.
  const CHORDS = { I: [0, 2, 4], ii: [1, 3, 5], IV: [3, 5, 0], V: [4, 6, 1], vi: [5, 0, 2] };
  const PROGRESSIONS = [['I', 'IV', 'V', 'I'], ['I', 'vi', 'IV', 'V'], ['I', 'ii', 'V', 'I'], ['I', 'V', 'vi', 'IV'], ['I', 'IV', 'I', 'V']];
  const mod7 = x => ((x % 7) + 7) % 7;
  const semisOf = (s, scale) => 12 * Math.floor(s / 7) + scale[mod7(s)];
  const isTritone = (a, b, scale) => Math.abs(semisOf(a, scale) - semisOf(b, scale)) % 12 === 6;

  function weightedPick(cands) {
    let r = Math.random() * cands.reduce((sum, [, w]) => sum + w, 0);
    for (const [value, w] of cands) if ((r -= w) <= 0) return value;
    return cands[cands.length - 1][0];
  }

  function scaleSteps(L, scale) {
    const n = L.notes;
    for (let tries = 0; tries < 500; tries++) {
      const prog = pick(PROGRESSIONS);
      const chordAt = i => CHORDS[prog[Math.min(prog.length - 1, Math.floor(i * prog.length / (n - 2)))]];
      const shape = pick(['arch', 'arch', 'valley', 'rise']);
      const contour = i => {
        const t = i / (n - 1);
        if (shape === 'arch') return L.arch * Math.sin(Math.PI * t);
        if (shape === 'valley') return -0.7 * L.arch * Math.sin(Math.PI * t);
        return 0.8 * L.arch * t;
      };

      const s = [pick([0, 2, 4, -3].filter(x => x >= L.low && x <= L.high))];
      let stuck = false;
      for (let i = 1; i < n - 2; i++) {
        const prev = s[i - 1];
        const last = i >= 2 ? prev - s[i - 2] : 0;
        const prevIsPassing = i >= 2 && !chordAt(i - 1).includes(mod7(prev));
        let run = 0;
        for (let j = i - 1; j >= 1 && last !== 0 && Math.sign(s[j] - s[j - 1]) === Math.sign(last); j--) run++;

        const cands = [];
        for (let c = L.low; c <= L.high; c++) {
          const iv = c - prev, size = Math.abs(iv);
          let w = L.leaps[size] || 0;
          if (!w || (iv === 0 && last === 0)) continue;                    // no three repeats in a row
          if (!L.tritone && isTritone(c, prev, L.scale)) continue;
          if (Math.abs(last) >= 3) {                                        // after a leap, step back
            w *= Math.sign(iv) === -Math.sign(last) && size <= 2 ? 4 : Math.sign(iv) === Math.sign(last) ? 0.15 : 1;
          }
          if (prevIsPassing && size !== 1) w *= 0.1;                        // passing notes resolve by step
          const chordTone = chordAt(i).includes(mod7(c));
          if (i % 2 === 0) w *= chordTone ? 3 : 0.3;                        // strong beat: chord tones
          else if (!chordTone && size !== 1) w *= 0.2;                      // weak beat: non-chord tones by step
          if (run >= 3 && Math.sign(iv) === Math.sign(last)) w *= 0.35;     // don't run on forever
          w *= Math.exp(-Math.abs(c - (s[0] + contour(i))) / 2.2);         // follow the contour
          cands.push([c, w]);
        }
        if (!cands.length) { stuck = true; break; }
        s.push(weightedPick(cands));
      }
      if (stuck) continue;

      // Cadence: the tonic nearest the line, approached by step from above or below
      const a = s[n - 3];
      const fin = [-7, 0, 7].filter(x => x >= L.low && x <= L.high)
        .reduce((best, x) => (Math.abs(x - a) < Math.abs(best - a) ? x : best));
      const pens = [fin + 1, fin - 1].filter(x => x >= L.low && x <= L.high && x !== a)
        .sort((x, y) => Math.abs(x - a) - Math.abs(y - a));
      if (!pens.length) continue;
      const pen = pens[0];
      if (Math.abs(pen - a) > 4 || (!L.tritone && isTritone(pen, a, L.scale))) continue;
      s.push(pen, fin);

      if (Math.max(...s) - Math.min(...s) > L.span) continue;
      if (new Set(s).size < Math.min(4, n - 1)) continue;                  // too static
      return s;
    }
    return null;
  }

  // Turn scale steps into spelled notes in a register that fits the staff and keyboard.
  function realize(key, L, s) {
    const scale = key.scale;
    const raise = s.map(() => 0);
    if (key.mode === 'minor') {
      for (let i = 0; i < s.length - 1; i++) {
        if (mod7(s[i]) === 6 && s[i + 1] === s[i] + 1) {                   // leading tone up to the tonic
          raise[i] = 1;
          if (i > 0 && s[i - 1] === s[i] - 1) raise[i - 1] = 1;           // melodic minor: no augmented 2nd
        }
      }
    }
    if (L.chromatic) {                                                     // chromatic lower neighbours: G F G -> G F♯ G
      for (let i = 1; i < s.length - 2; i++) {
        if (s[i - 1] === s[i + 1] && s[i] === s[i - 1] - 1 && !raise[i]
          && semisOf(s[i - 1], scale) - semisOf(s[i], scale) === 2 && Math.random() < L.chromatic) raise[i] = 1;
      }
    }
    const rel = s.map((x, i) => semisOf(x, scale) + raise[i]);
    const registers = [];
    for (let oct = 2; oct <= 5; oct++) {
      const td = dia(key.letter, oct), tm = diaMidi(td) + key.alter;
      const lo = tm + Math.min(...rel), hi = tm + Math.max(...rel);
      if (lo >= LOW && hi <= HIGH) registers.push({ td, tm, off: Math.abs((lo + hi) / 2 - 69) });
    }
    if (!registers.length) return null;
    const { td, tm } = registers.sort((a, b) => a.off - b.off)[0];
    const notes = s.map((x, i) => {
      const d = td + x, midi = tm + rel[i];
      return { d, midi, alter: midi - diaMidi(d) };
    });
    if (notes.some(x => Math.abs(x.alter) > 1)) return null;
    return { key, notes, tonicMidi: tm };
  }

  function makeMelody(L) {
    const keys = KEYS.filter(L.keys);
    for (;;) {
      const key = pick(keys);
      const s = scaleSteps({ ...L, scale: key.scale }, key.scale);
      const m = s && realize(key, L, s);
      if (m) return m;
    }
  }

  // Spell what the player types the way the key would: scale notes first,
  // then plain letters, then sharps (sharp keys) or flats (flat keys).
  function spellerFor(key) {
    const tl = LETTERS.indexOf(key.letter);
    const tonicPc = (SEMIS[tl] + key.alter + 12) % 12;
    const known = new Map();
    const add = (k, extra) => {
      const L = (tl + k) % 7;
      const pc = (tonicPc + key.scale[k] + extra + 12) % 12;
      const alter = ((pc - SEMIS[L] + 18) % 12) - 6;
      if (!known.has(pc) && Math.abs(alter) <= 1) known.set(pc, [L, alter]);
    };
    for (let k = 0; k < 7; k++) add(k, 0);
    if (key.mode === 'minor') { add(6, 1); add(5, 1); }
    return midi => {
      const pc = midi % 12;
      let spelled = known.get(pc);
      if (!spelled) {
        const nat = SEMIS.indexOf(pc);
        spelled = nat >= 0 ? [nat, 0] : key.sig < 0 ? [SEMIS.indexOf((pc + 1) % 12), -1] : [SEMIS.indexOf((pc + 11) % 12), 1];
      }
      const [L, alter] = spelled;
      return { d: (Math.floor((midi - alter) / 12) - 1) * 7 + L, midi, alter };
    };
  }

  // No key signature, so the melody reads as one long bar: an accidental lasts
  // until the same note changes again, which then gets its own sign.
  function withAccidentals(notes) {
    const inEffect = new Map();
    return notes.map(n => {
      const cur = inEffect.get(n.d) ?? 0;
      inEffect.set(n.d, n.alter);
      return { ...n, acc: n.alter !== cur ? accFor(n.alter) || 'natural' : '' };
    });
  }

  // ---------- Setup ----------
  const levelGrid = $('levelGrid'), countGroup = $('mCountGroup'), chordToggle = $('keyChordToggle');

  function buildSetup() {
    levelGrid.innerHTML = Object.entries(LEVELS).map(([id, L]) => `
      <button type="button" class="level-tile" role="radio" id="level-${id}" data-level="${id}" aria-checked="false">
        <span class="level-name">${L.name}</span>
        <span class="level-notes">${L.notes} notes</span>
        <span class="level-desc">${L.desc}</span>
      </button>`).join('');
    countGroup.innerHTML = COUNTS.map(n =>
      `<button type="button" role="radio" id="mcount-${n}" data-count="${n}" aria-checked="false">${n}</button>`).join('');

    levelGrid.addEventListener('click', e => {
      const b = e.target.closest('[data-level]');
      if (!b) return;
      settings.level = b.dataset.level;
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
    chordToggle.addEventListener('click', () => {
      settings.keyChord = !settings.keyChord;
      save();
      syncSetup();
    });
  }

  function syncSetup() {
    for (const b of levelGrid.querySelectorAll('[data-level]')) b.setAttribute('aria-checked', String(b.dataset.level === settings.level));
    for (const b of countGroup.querySelectorAll('[data-count]')) b.setAttribute('aria-checked', String(Number(b.dataset.count) === settings.count));
    chordToggle.setAttribute('aria-checked', String(settings.keyChord));
  }

  // ---------- Game ----------
  const stage = $('mStage'), card = $('mCard'), scroll = $('mScroll'), svg = $('mStaff');
  const caption = $('mCaption'), confirmBtn = $('mConfirm');
  const timers = App.makeTimers();
  let keyEls = [];
  let game = null;

  function playMelody(withChord) {
    const { m, level } = game;
    App.stopAll();   // replays never overlap
    let t = 0;
    if (withChord && settings.keyChord) {
      const root = m.tonicMidi >= 60 ? m.tonicMidi - 12 : m.tonicMidi;
      for (const k of [0, 2, 4]) App.play(root + m.key.scale[k], 0, 0.55, 1.1);
      t = 1.7;
    }
    m.notes.forEach((x, i) => App.play(x.midi, t + i * level.beat, 1, level.beat * 0.95));
  }

  async function start() {
    const loading = App.ensureAudio();
    timers.clear();
    game = { level: LEVELS[settings.level], i: 0, total: settings.count, correct: 0, notes: 0, perfect: 0, times: [], m: null, phase: 'loading' };
    App.show('melodyGame');
    caption.textContent = 'Loading piano…';
    confirmBtn.disabled = true;
    App.setProgress($('mProgressFill'), $('mProgressText'), 0, game.total);
    const g = game;
    await Promise.race([loading, App.wait(5000)]);
    if (game !== g) return;   // left the screen while loading
    next();
  }

  function next() {
    const m = makeMelody(game.level);
    Object.assign(game, {
      m,
      spell: spellerFor(m.key),
      input: [m.notes[0]],     // the first note is given
      phase: 'input',
      octave: Math.min(OCT_MAX, Math.max(OCT_MIN, Math.floor(m.notes[0].midi / 12) - 1)),
      start: 0,
    });
    caption.textContent = `${m.key.name} · ${m.notes.length} notes · first note given`;
    App.setProgress($('mProgressFill'), $('mProgressText'), game.i, game.total);
    update();
    draw(true);
    timers.later(() => {
      if (!game) return;
      playMelody(true);
      game.start = performance.now();
    }, 300);
  }

  function update() {
    const n = game.m ? game.m.notes.length : 0;
    $('mOctave').textContent = `C${game.octave}–B${game.octave}`;
    $('mDown').disabled = game.octave <= OCT_MIN;
    $('mUp').disabled = game.octave >= OCT_MAX;
    $('mCount').textContent = game.m ? `${game.input.length} / ${n}` : '';
    $('mDelete').disabled = game.phase !== 'input' || game.input.length <= 1;
    if (game.phase === 'review') {
      confirmBtn.textContent = game.i + 1 >= game.total ? 'See results' : 'Next melody';
      confirmBtn.disabled = false;
    } else {
      confirmBtn.textContent = 'Confirm';
      confirmBtn.disabled = game.phase !== 'input' || game.input.length < n;
    }
  }

  function draw(instant) {
    if (!game || !game.m) return;
    const { m, input, phase } = game;
    const n = m.notes.length;
    const mine = withAccidentals(input), truth = withAccidentals(m.notes);
    const cols = [];
    for (let i = 0; i < n; i++) {
      if (phase !== 'review') {
        cols.push(i < input.length ? { notes: [{ ...mine[i], state: '' }] } : { notes: [], cursor: i === input.length });
      } else if (input[i].midi === m.notes[i].midi) {
        cols.push({ notes: [{ ...mine[i], state: i ? 'ok' : '' }] });
      } else {
        cols.push({ notes: [{ ...mine[i], state: 'bad' }, { ...truth[i], state: 'ok' }] });
      }
    }
    const layout = App.fitInStage(stage, card, scroll, maxHeight =>
      App.renderGrandStaff(svg, cols, { trebleOnly: true, cols: n, colWidth: 60, fit: scroll.clientWidth, min: 0.8, max: 1.5, maxHeight }));
    const focus = phase === 'review' ? 0 : Math.min(input.length, n - 1);
    scroll.scrollTo({
      left: Math.max(0, layout.centers[focus] * layout.scale - scroll.clientWidth * 0.55),
      behavior: instant || reduceMotion.matches ? 'auto' : 'smooth',
    });
  }

  function press(pc) {
    App.ensureAudio();
    if (!game || App.current() !== 'melodyGame') return;
    const midi = 12 * (game.octave + 1) + pc;
    App.play(midi, 0, 1, 0.7);
    if (game.phase !== 'input' || game.input.length >= game.m.notes.length) return;
    game.input.push(game.spell(midi));
    update();
    draw();
  }

  function removeLast() {
    if (!game || game.phase !== 'input' || game.input.length <= 1) return;
    game.input.pop();
    update();
    draw();
  }

  function shiftOctave(by) {
    if (!game) return;
    game.octave = Math.min(OCT_MAX, Math.max(OCT_MIN, game.octave + by));
    update();
  }

  function confirm() {
    if (!game || !game.m) return;
    if (game.phase === 'review') return advance();
    const { m, input } = game;
    const n = m.notes.length;
    if (game.phase !== 'input' || input.length < n) return;
    game.phase = 'review';
    game.times.push(performance.now() - (game.start || performance.now()));
    // The given first note doesn't count
    const right = input.filter((x, i) => i > 0 && x.midi === m.notes[i].midi).length;
    game.correct += right;
    game.notes += n - 1;
    if (right === n - 1) game.perfect++;
    caption.textContent = right === n - 1 ? 'Every note correct' : `${right} of ${n - 1} notes correct · the right notes are in green`;
    update();
    draw();
    if (right < n - 1) timers.later(() => game && playMelody(false), 400);   // hear the real one again
  }

  function advance() {
    game.i++;
    if (game.i >= game.total) return finish();
    next();
  }

  function finish() {
    const g = game;
    App.showResults({
      correct: g.correct,
      total: g.notes,
      times: g.times,
      unit: 'melody',
      statsKey: 'melody',
      summary: `${g.correct} of ${g.notes} notes correct · ${g.perfect} of ${g.total} ${g.total === 1 ? 'melody' : 'melodies'} perfect.`,
      review: null,
      onAgain: start,
      onSettings: openSetup,
    });
    game = null;
  }

  function openSetup() {
    timers.clear();
    game = null;
    syncSetup();
    App.show('melodySetup');
  }

  // ---------- Wiring ----------
  buildSetup();
  syncSetup();
  keyEls = App.makePiano($('mKeys'), press);
  $('mStartBtn').addEventListener('click', start);
  $('mQuitBtn').addEventListener('click', openSetup);
  $('mReplay').addEventListener('click', () => { App.ensureAudio(); if (game && game.m) playMelody(false); });
  $('mDown').addEventListener('click', () => shiftOctave(-1));
  $('mUp').addEventListener('click', () => shiftOctave(1));
  $('mDelete').addEventListener('click', removeLast);
  confirmBtn.addEventListener('click', confirm);

  // Computer keyboard: A S D F G H J = white keys, W E T Y U = black keys
  const KEYMAP = { a: 0, w: 1, s: 2, e: 3, d: 4, f: 5, t: 6, g: 7, y: 8, h: 9, u: 10, j: 11 };
  document.addEventListener('keydown', e => {
    if (e.repeat || e.metaKey || e.ctrlKey || e.altKey) return;
    if (App.current() !== 'melodyGame') return;
    const k = e.key.toLowerCase();
    if (k === 'escape') return openSetup();
    if (k === 'enter') { e.preventDefault(); return confirm(); }
    if (k === 'backspace') { e.preventDefault(); return removeLast(); }
    if (k === ' ') { e.preventDefault(); if (game && game.m) playMelody(false); return; }
    if (k === 'z') return shiftOctave(-1);
    if (k === 'x') return shiftOctave(1);
    if (k in KEYMAP) {
      e.preventDefault();
      App.tapKey(keyEls[KEYMAP[k]]);
      press(KEYMAP[k]);
    }
  });

  App.watchStage(stage, () => { if (App.current() === 'melodyGame') draw(true); });
  App.onRedraw(() => { if (App.current() === 'melodyGame') draw(true); });

  window.Melody = { open: openSetup };
})();
