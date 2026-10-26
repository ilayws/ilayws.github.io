// Perfect Pitch drill: hear one note with nothing to compare it to, and find it on a one-octave piano.
(() => {
  'use strict';
  const { CLEFS, $, randInt } = App;

  // How each pitch class is named and written: [letter index, alteration]
  const SPELL = [[0, 0], [0, 1], [1, 0], [2, -1], [2, 0], [3, 0], [3, 1], [4, 0], [5, -1], [5, 0], [6, -1], [6, 0]];
  const SHORT = ['C', 'C♯', 'D', 'E♭', 'E', 'F', 'F♯', 'G', 'A♭', 'A', 'B♭', 'B'];
  const FULL = ['C', 'C♯ / D♭', 'D', 'D♯ / E♭', 'E', 'F', 'F♯ / G♭', 'G', 'G♯ / A♭', 'A', 'A♯ / B♭', 'B'];
  const LOW = 48, HIGH = 83;             // C3 to B5
  const COUNTS = [5, 10, 20, 30];
  const RIGHT_PAUSE = 900, RETRY_DELAY = 700;
  const KEYMAP = { a: 0, w: 1, s: 2, e: 3, d: 4, f: 5, t: 6, g: 7, y: 8, h: 9, u: 10, j: 11 };
  const KEY = 'pitch.settings';

  const settings = App.loadSettings(KEY, { count: 10 }, s => COUNTS.includes(s.count));
  settings.endless = settings.endless === true;
  const save = () => App.saveSettings(KEY, settings);

  const card = $('ppCard'), prompt = $('ppPrompt'), nameEl = $('ppName'), sub = $('ppSub');
  const replayBtn = $('ppReplay'), nextBtn = $('ppNext');
  const fill = $('ppProgressFill'), text = $('ppProgressText');
  const timers = App.makeTimers();

  // ---------- Setup ----------
  const syncCount = App.countPicker($('ppCountGroup'), COUNTS, settings, () => { save(); syncCount(); }, 'ppcount');

  // ---------- Notes ----------
  // d: diatonic number, alter: sharp or flat, as the note is written
  function spell(midi) {
    const [letter, alter] = SPELL[midi % 12];
    return { d: (Math.floor(midi / 12) - 1) * 7 + letter, alter };
  }
  const written = midi => { const { d, alter } = spell(midi); return App.noteName(d, alter); };

  // A new pitch class each time, so two notes in a row never match
  function makeNote(prev) {
    let midi;
    do midi = randInt(LOW, HIGH); while (prev && midi % 12 === prev.pc);
    return { midi, pc: midi % 12, start: 0 };
  }

  // The key that was pressed, in the octave nearest the note
  const nearest = (pc, midi) => midi + ((pc - (midi % 12) + 18) % 12) - 6;

  // ---------- Piano ----------
  const keyEls = App.makePiano($('ppKeys'), pc => press(pc));
  const clearKeys = () => keyEls.forEach(el => el.classList.remove('ok', 'bad', 'down'));

  // ---------- Game ----------
  // state: 'loading' | 'listen' (waiting for an answer) | 'right' | 'missed' (waiting for Continue)
  let game = null;

  function hear(midi) {
    App.stopAll();
    App.play(midi);
  }

  function showProgress() {
    const right = game.results.filter(r => r === 'ok').length;
    if (game.endless) App.setRun(fill, text, game.lives, right);
    else App.setProgress(fill, text, game.i, game.count);
  }

  function setCard(state, promptText, name = '', subText = '') {
    card.dataset.state = state;
    prompt.textContent = promptText;
    nameEl.textContent = name;
    sub.textContent = subText;
  }

  async function start() {
    const loading = App.ensureAudio();
    timers.clear();
    clearKeys();
    game = {
      endless: settings.endless,
      count: settings.count,
      lives: App.LIVES,
      i: 0,
      q: null,
      notes: [],
      results: [],
      times: [],
      misses: [],
      state: 'loading',
    };
    nextBtn.hidden = true;
    App.show('pitchGame');
    showProgress();
    setCard('listen', 'Loading piano…');
    const g = game;
    await Promise.race([loading, App.wait(5000)]);
    if (game !== g) return;   // left the screen while loading
    nextNote();
  }

  function nextNote() {
    clearKeys();
    nextBtn.hidden = true;
    game.q = makeNote(game.q);
    game.state = 'listen';
    showProgress();
    setCard('listen', 'Which note is this?');
    timers.later(() => {
      if (!game || game.state !== 'listen') return;
      hear(game.q.midi);
      game.q.start = performance.now();
    }, 300);
  }

  function press(pc) {
    App.ensureAudio();
    if (App.current() !== 'pitchGame' || !game || !game.q) return;
    const q = game.q;
    if (game.state === 'missed') return hear(nearest(pc, q.midi));   // comparing keys after a miss
    if (game.state !== 'listen' || !q.start) return;

    const ms = performance.now() - q.start;
    const ok = pc === q.pc;
    game.notes.push(q.midi);
    game.results.push(ok ? 'ok' : 'bad');
    game.times.push(ms);
    App.recordItem('pitch', String(q.pc), ok, ms, { label: SHORT[q.pc] });
    hear(nearest(pc, q.midi));

    if (ok) {
      game.state = 'right';
      keyEls[pc].classList.add('ok');
      setCard('ok', 'Right', FULL[q.pc], written(q.midi));
      showProgress();
      timers.later(advance, RIGHT_PAUSE);
      return;
    }
    App.recordMixup('pitch', String(q.pc), String(pc));
    game.misses.push([q.pc, pc]);
    game.lives--;
    game.state = 'missed';
    keyEls[pc].classList.add('bad');
    keyEls[q.pc].classList.add('ok');
    setCard('bad', 'It was', FULL[q.pc], `${written(q.midi)} · you played ${FULL[pc]}`);
    showProgress();
    timers.later(() => { if (game && game.q === q) hear(q.midi); }, RETRY_DELAY);
    nextBtn.textContent = isLast() ? 'See results' : 'Continue';
    nextBtn.hidden = false;
  }

  const isLast = () => (game.endless ? game.lives <= 0 : game.i + 1 >= game.count);

  function advance() {
    if (!game) return;
    if (isLast()) return finish();
    game.i++;
    nextNote();
  }

  function continueAfterMiss() {
    if (!game || game.state !== 'missed') return;
    timers.clear();
    advance();
  }

  function replay() {
    App.ensureAudio();
    if (game && game.q && game.q.start && game.state !== 'right') hear(game.q.midi);
  }

  // Misses a half step away are near misses; say how many there were
  function describeMisses(misses) {
    if (!misses.length) return '';
    const near = misses.filter(([a, g]) => [1, 11].includes((a - g + 12) % 12)).length;
    if (!near) return 'None of your misses were a half step off.';
    return `${near} of ${misses.length} ${misses.length === 1 ? 'miss was' : 'misses were'} a half step off.`;
  }

  // The notes you heard, on treble and bass staves; one long measure, so accidentals carry
  function reviewItems(notes, results) {
    const items = [], carried = new Map();
    let prevClef = null;
    notes.forEach((midi, i) => {
      const clef = midi >= 60 ? 'treble' : 'bass';
      if (clef !== prevClef) { items.push({ type: 'clef', clef }); prevClef = clef; }
      const { d, alter } = spell(midi);
      const inEffect = carried.get(d) ?? 0;
      const acc = alter !== inEffect ? App.accFor(alter) || 'natural' : null;
      carried.set(d, alter);
      items.push({ type: 'note', p: d - CLEFS[clef].bottomLine, acc, index: i, state: results[i] });
    });
    return items;
  }

  function finish() {
    const done = game;
    timers.clear();
    game = null;
    App.showResults({
      correct: done.results.filter(r => r === 'ok').length,
      total: done.results.length,
      times: done.times,
      unit: 'note',
      statsKey: 'pitch',
      endless: done.endless ? { level: 'any' } : null,
      extra: describeMisses(done.misses),
      review: (svg, fit) => App.renderStaff(svg, reviewItems(done.notes, done.results), { fit, min: 0.6, max: 1.2 }),
      onAgain: start,
      onSettings: openSetup,
    });
  }

  function openSetup() {
    timers.clear();
    clearKeys();
    game = null;
    syncCount();
    App.show('pitchSetup');
  }

  // ---------- Wiring ----------
  syncCount();
  $('ppStartBtn').addEventListener('click', start);
  $('ppQuitBtn').addEventListener('click', openSetup);
  replayBtn.addEventListener('click', replay);
  nextBtn.addEventListener('click', continueAfterMiss);

  document.addEventListener('keydown', e => {
    if (e.repeat || e.metaKey || e.ctrlKey || e.altKey) return;
    if (App.current() !== 'pitchGame') return;
    const k = e.key.toLowerCase();
    if (k === 'escape') return openSetup();
    if (k === 'enter') { e.preventDefault(); return continueAfterMiss(); }
    if (k === ' ') { e.preventDefault(); return replay(); }
    if (k in KEYMAP) {
      e.preventDefault();
      App.tapKey(keyEls[KEYMAP[k]]);
      press(KEYMAP[k]);
    }
  });

  window.Pitch = {
    open: openSetup,
    // For the stats page: the twelve notes from C, written above middle C
    order: SHORT.map((_, pc) => String(pc)),
    label: id => SHORT[Number(id)],
    shape: id => [{ ...spell(60 + Number(id)), midi: 60 + Number(id) }],
  };
})();
