// Clefs drill: read notes across treble, bass and alto clef and play them on a one-octave piano.
(() => {
  'use strict';
  const { LETTERS, SEMIS, CLEFS, diaName, diaMidi, $, randInt, pick, reduceMotion } = App;

  const P_MIN = -6, P_MAX = 14;          // up to three ledger lines below / above
  const ACC_CHANCE = 0.4;
  // Skip E#, B#, Fb and Cb: they land on white keys and mostly confuse.
  const NO_SHARP = new Set([2, 6]), NO_FLAT = new Set([0, 3]);
  const COUNTS = [5, 10, 20, 30];
  const WRONG_PAUSE = 1000;
  const LOOKAHEAD = 10;                  // endless: notes kept on the staff ahead of the current one
  const KEY = 'clefs.settings';

  const settings = App.loadSettings(KEY, { clefs: ['treble', 'bass'], count: 10, accidentals: false },
    s => Array.isArray(s.clefs) && COUNTS.includes(s.count));
  settings.clefs = settings.clefs.filter(c => c in CLEFS);
  settings.endless = settings.endless === true;
  const save = () => App.saveSettings(KEY, settings);

  const clefGrid = $('clefGrid'), countGroup = $('countGroup');
  const startBtn = $('startBtn'), hint = $('hint'), accToggle = $('accToggle');
  const staffSvg = $('staff'), staffScroll = $('staffScroll');
  const keysEl = $('keys');
  const stageEl = $('clefsStage'), cardEl = $('clefsCard');
  const timers = App.makeTimers();

  // ---------- Setup ----------
  const syncCount = App.countPicker(countGroup, COUNTS, settings, () => { save(); syncSetup(); }, 'count');

  function buildSetup() {
    clefGrid.innerHTML = Object.entries(CLEFS).map(([k, c]) => `
      <button type="button" class="clef-tile" id="clef-${k}" data-clef="${k}" aria-pressed="false">
        <svg class="mini" aria-hidden="true"></svg>
        <span class="name">${c.name}</span>
        <span class="range">${diaName(c.bottomLine + P_MIN)}–${diaName(c.bottomLine + P_MAX)}</span>
      </button>`).join('');

    clefGrid.addEventListener('click', e => {
      const tile = e.target.closest('[data-clef]');
      if (!tile) return;
      const k = tile.dataset.clef;
      settings.clefs = settings.clefs.includes(k)
        ? settings.clefs.filter(c => c !== k)
        : Object.keys(CLEFS).filter(c => c === k || settings.clefs.includes(c));
      save();
      syncSetup();
    });
    accToggle.addEventListener('click', () => {
      settings.accidentals = !settings.accidentals;
      save();
      syncSetup();
    });
  }

  function syncSetup() {
    for (const tile of clefGrid.querySelectorAll('[data-clef]')) {
      tile.setAttribute('aria-pressed', String(settings.clefs.includes(tile.dataset.clef)));
    }
    syncCount();
    accToggle.setAttribute('aria-checked', String(settings.accidentals));
    const none = settings.clefs.length === 0;
    startBtn.disabled = none;
    hint.hidden = !none;
  }

  function renderSetupStaves() {
    for (const tile of clefGrid.querySelectorAll('[data-clef]')) {
      App.renderStaff(tile.querySelector('svg'), [{ type: 'clef', clef: tile.dataset.clef }], { scale: 0.62, barline: false, range: [11.6, -4] });
    }
  }

  // ---------- Piano ----------
  let keyEls = [];
  const buildPiano = () => { keyEls = App.makePiano(keysEl, pc => press(pc)); };

  // ---------- Game ----------
  // Runs of 1–4 notes per clef; the clef changes between runs when more than one is picked.
  // The round reads as one long measure, so an accidental carries forward to later
  // notes of the same pitch; those get a natural sign when they go back to natural.
  // noteStream returns a function that gives the next note each time it's called.
  function noteStream(clefs, withAcc) {
    const carried = new Map();   // diatonic note -> alteration still in effect
    let clef = pick(clefs), run = randInt(1, 4), prev = null;
    return () => {
      if (run === 0) {
        if (clefs.length > 1) clef = pick(clefs.filter(c => c !== clef));
        run = randInt(1, 4);
      }
      run--;
      const base = CLEFS[clef].bottomLine;
      let p, d, alter, midi;
      do {
        p = randInt(P_MIN, P_MAX);
        d = base + p;
        alter = 0;
        if (withAcc && Math.random() < ACC_CHANCE) {
          const opts = [];
          if (!NO_SHARP.has(d % 7)) opts.push(1);
          if (!NO_FLAT.has(d % 7)) opts.push(-1);
          alter = pick(opts);
        }
        midi = diaMidi(d) + alter;
      } while (midi === prev);
      prev = midi;

      const inEffect = carried.get(d) ?? 0;
      let acc = null;
      if (alter === 1) acc = 'sharp';
      else if (alter === -1) acc = 'flat';
      else if (inEffect !== 0) acc = 'natural';
      carried.set(d, alter);
      return { clef, p, acc, midi };
    };
  }

  function itemsFor(notes, results, current) {
    const items = [];
    let prevClef = null;
    notes.forEach((n, i) => {
      if (n.clef !== prevClef) { items.push({ type: 'clef', clef: n.clef }); prevClef = n.clef; }
      items.push({ type: 'note', p: n.p, acc: n.acc, index: i, state: results[i] || (i === current ? 'current' : '') });
    });
    return items;
  }

  let game = null;

  function clearFeedback() {
    timers.clear();
    keyEls.forEach(el => el.classList.remove('ok', 'bad', 'down'));
  }

  function start() {
    App.ensureAudio();
    clearFeedback();
    const next = noteStream(settings.clefs, settings.accidentals);
    const endless = settings.endless;
    game = {
      endless,
      next,
      lives: App.LIVES,
      notes: Array.from({ length: endless ? LOOKAHEAD + 1 : settings.count }, next),
      i: 0,
      results: [],
      times: [],
      locked: false,
      noteStart: 0,
    };
    App.show('clefsGame');
    drawGame(true);
    game.noteStart = performance.now();
  }

  function drawGame(instant) {
    const { notes, i } = game;
    const items = itemsFor(notes, game.results, i);
    const layout = App.fitInStage(stageEl, cardEl, staffScroll, maxHeight =>
      App.renderStaff(staffSvg, items, { fit: staffScroll.clientWidth, maxHeight, min: 1, max: 1.6 }));
    if (game.endless) App.setRun($('progressFill'), $('progressText'), game.lives, game.results.filter(r => r === 'ok').length);
    else App.setProgress($('progressFill'), $('progressText'), i, notes.length);

    const x = layout.noteX[Math.min(i, notes.length - 1)] * layout.scale;
    staffScroll.scrollTo({
      left: Math.max(0, x - staffScroll.clientWidth * 0.4),
      behavior: instant || reduceMotion.matches ? 'auto' : 'smooth',
    });
  }

  function press(pc) {
    App.ensureAudio();
    if (App.current() !== 'clefsGame' || !game || game.locked) return;

    const note = game.notes[game.i];
    const targetPc = note.midi % 12;
    App.play(note.midi + ((pc - targetPc + 18) % 12) - 6);   // the pressed key nearest the written note
    const ms = performance.now() - game.noteStart;
    game.times.push(ms);
    const d = CLEFS[note.clef].bottomLine + note.p, alter = note.midi - diaMidi(d);
    const name = App.noteName(d, alter);
    App.recordItem('clefs', `${note.clef}:${name}`, pc === targetPc, ms, { clef: note.clef, p: note.p, label: name });

    if (pc === targetPc) {
      game.results[game.i] = 'ok';
      keyEls[pc].classList.add('ok');
      timers.later(() => keyEls[pc].classList.remove('ok'), 220);
      advance();
    } else {
      game.results[game.i] = 'bad';
      game.locked = true;
      game.lives--;
      keyEls[pc].classList.add('bad');
      keyEls[targetPc].classList.add('ok');
      drawGame();
      timers.later(() => {
        keyEls[pc].classList.remove('bad');
        keyEls[targetPc].classList.remove('ok');
        if (game.endless && game.lives <= 0) return finish();
        game.locked = false;
        advance();
      }, WRONG_PAUSE);
    }
  }

  function advance() {
    game.i++;
    if (game.endless) while (game.notes.length - game.i <= LOOKAHEAD) game.notes.push(game.next());
    drawGame();
    if (game.i >= game.notes.length) {
      game.locked = true;
      timers.later(finish, 450);
    } else {
      game.noteStart = performance.now();
    }
  }

  function finish() {
    const done = game;
    done.notes = done.notes.slice(0, done.results.length);   // endless: drop the notes still ahead
    App.showResults({
      correct: done.results.filter(r => r === 'ok').length,
      total: done.notes.length,
      times: done.times,
      unit: 'note',
      statsKey: 'clefs',
      endless: done.endless ? { level: 'any' } : null,
      review: (svg, fit) => App.renderStaff(svg, itemsFor(done.notes, done.results, -1), { fit, min: 0.6, max: 1.2 }),
      onAgain: start,
      onSettings: openSetup,
    });
  }

  function openSetup() {
    clearFeedback();
    game = null;
    syncSetup();
    App.show('clefsSetup');
  }

  // ---------- Wiring ----------
  buildSetup();
  buildPiano();
  syncSetup();
  startBtn.addEventListener('click', start);
  $('quitBtn').addEventListener('click', openSetup);

  document.addEventListener('keydown', e => {
    if (e.repeat || e.metaKey || e.ctrlKey || e.altKey) return;
    if (App.current() !== 'clefsGame') return;
    if (e.key === 'Escape') return openSetup();
    const idx = LETTERS.indexOf(e.key.toUpperCase());
    if (idx < 0) return;
    e.preventDefault();
    App.tapKey(keyEls[SEMIS[idx]]);
    press(SEMIS[idx]);
  });

  App.watchStage(stageEl, () => { if (App.current() === 'clefsGame' && game) drawGame(true); });
  App.onRedraw(() => {
    renderSetupStaves();
    if (App.current() === 'clefsGame' && game) drawGame(true);
  });

  window.Clefs = { open: openSetup };
})();
