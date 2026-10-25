// Progressions drill: hear a chord progression and name each chord as a scale degree.
(() => {
  'use strict';
  const { LETTERS, SEMIS, $, pick, reduceMotion } = App;

  // ---------- Chord symbols ----------
  // Numerals are relative to the major scale on the tonic, so in a minor key the
  // chord on the 6th note is ♭VI. Case shows quality: upper = major, lower = minor.
  // Suffixes: ° diminished, + augmented, 7 (dominant on upper case, minor 7th on
  // lower case), maj7, ø7 half-diminished, °7 diminished 7th.
  const MAJOR = [0, 2, 4, 5, 7, 9, 11];
  const ROMANS = ['I', 'II', 'III', 'IV', 'V', 'VI', 'VII'];
  const QUALITIES = {
    '': up => (up ? [0, 4, 7] : [0, 3, 7]),
    '°': () => [0, 3, 6],
    '+': () => [0, 4, 8],
    '7': up => (up ? [0, 4, 7, 10] : [0, 3, 7, 10]),
    maj7: () => [0, 4, 7, 11],
    'ø7': () => [0, 3, 6, 10],
    '°7': () => [0, 3, 6, 9],
  };
  const TONE_STEPS = [0, 2, 4, 6];          // root, 3rd, 5th, 7th in letter names

  function parse(sym) {
    const m = /^([b#]?)(VII|VI|V|IV|III|II|I|vii|vi|v|iv|iii|ii|i)(°7|ø7|°|\+|maj7|7)?$/.exec(sym);
    if (!m) throw new Error(`Unknown chord symbol ${sym}`);
    const [, acc, roman, q = ''] = m;
    const up = roman === roman.toUpperCase();
    return {
      sym,
      degree: ROMANS.indexOf(roman.toUpperCase()),
      shift: acc === 'b' ? -1 : acc === '#' ? 1 : 0,
      semis: QUALITIES[q](up),
      // What the player answers: degree and triad quality, without the 7th
      base: acc + roman + (q === '°' || q === 'ø7' || q === '°7' ? '°' : q === '+' ? '+' : ''),
    };
  }
  const pretty = s => s.replace(/^b/, '♭').replace(/^#/, '♯');

  // ---------- Library ----------
  const P = (name, style, chords) => ({ name, style, chords });
  const LEVELS = {
    easy: {
      name: 'Easy', desc: 'Major keys · plain triads · 3–5 chords', beat: 1.25,
      pad: ['I', 'ii', 'iii', 'IV', 'V', 'vi'],
      list: [
        P('Authentic cadence', 'Classical', ['I', 'IV', 'V', 'I']),
        P('Tonic and dominant', 'Folk', ['I', 'V', 'I']),
        P('Plagal "Amen" cadence', 'Classical, hymns', ['I', 'IV', 'I']),
        P('IV–V–I', 'Classical', ['IV', 'V', 'I']),
        P('ii–V–I in triads', 'Classical, jazz', ['I', 'ii', 'V', 'I']),
        P('Four-chord song', 'Pop', ['I', 'V', 'vi', 'IV']),
        P('50s doo-wop', 'Pop', ['I', 'vi', 'IV', 'V']),
        P('Pop, starting on vi', 'Pop', ['vi', 'IV', 'I', 'V']),
        P('I–IV–vi–V', 'Pop', ['I', 'IV', 'vi', 'V']),
        P('Turnaround in triads', 'Jazz, pop', ['I', 'vi', 'ii', 'V']),
        P('"La Bamba" riff', 'Rock', ['I', 'IV', 'V', 'IV']),
        P('I–iii–IV–V', 'Pop', ['I', 'iii', 'IV', 'V']),
        P('Descending to the cadence', 'Classical', ['I', 'vi', 'IV', 'V', 'I']),
        P('ii–V–I', 'Classical, jazz', ['ii', 'V', 'I']),
        P('Rock plagal cadence', 'Rock', ['V', 'IV', 'I']),
        P('I–V–IV–V', 'Rock, country', ['I', 'V', 'IV', 'V']),
        P('I–ii–IV–V', 'Pop', ['I', 'ii', 'IV', 'V']),
        P('I–IV–ii–V', 'Pop, jazz', ['I', 'IV', 'ii', 'V']),
        P('Circle into the tonic', 'Classical', ['vi', 'ii', 'V', 'I']),
        P('I–iii–vi–IV', 'Pop', ['I', 'iii', 'vi', 'IV']),
        P('Starting on IV', 'Pop', ['IV', 'I', 'V', 'vi']),
        P('vi–V–IV–V', 'Pop, rock', ['vi', 'V', 'IV', 'V']),
        P('Royal road in triads', 'J-pop', ['IV', 'V', 'iii', 'vi']),
        P('Hymn cadence', 'Hymns', ['I', 'IV', 'I', 'V', 'I']),
        P('Canon opening', 'Classical', ['I', 'V', 'vi', 'iii', 'IV']),
      ],
    },
    medium: {
      name: 'Medium', desc: 'Sevenths, minor keys, borrowed chords', beat: 1.1,
      pad: ['I', 'i', 'ii', 'II', 'bIII', 'iii', 'IV', 'iv', 'V', 'v', 'bVI', 'vi', 'bVII', 'vii°'],
      list: [
        P('ii–V–I', 'Jazz', ['ii7', 'V7', 'Imaj7']),
        P('I–vi–ii–V turnaround', 'Jazz', ['Imaj7', 'vi7', 'ii7', 'V7']),
        P('Minor cadence', 'Classical', ['i', 'iv', 'V', 'i']),
        P('Natural-minor cadence', 'Folk', ['i', 'iv', 'v', 'i']),
        P('Andalusian cadence', 'Classical, flamenco', ['i', 'bVII', 'bVI', 'V']),
        P('Epic minor loop', 'Pop', ['i', 'bVI', 'bIII', 'bVII']),
        P('Aeolian vamp', 'Rock', ['i', 'bVII', 'bVI', 'bVII']),
        P('Mixolydian rock', 'Rock', ['I', 'bVII', 'IV', 'I']),
        P('Minor plagal cadence', 'Pop', ['I', 'IV', 'iv', 'I']),
        P('Aeolian cadence', 'Pop, games', ['I', 'bVI', 'bVII', 'I']),
        P('"Eight Days a Week"', 'Pop', ['I', 'II7', 'IV', 'I']),
        P('Royal road', 'J-pop', ['IVmaj7', 'V7', 'iii7', 'vi']),
        P('Short blues', 'Blues, jazz', ['I7', 'IV7', 'I7', 'V7', 'IV7', 'I7']),
        P('Cadence with a dominant 7th', 'Classical', ['I', 'vi', 'IV', 'V7', 'I']),
        P("Pachelbel's Canon", 'Classical', ['I', 'V', 'vi', 'iii', 'IV', 'I', 'IV', 'V']),
        P('ii–V–I–vi', 'Jazz', ['ii7', 'V7', 'Imaj7', 'vi7']),
        P('ii–V–iii–vi', 'Jazz, pop', ['ii7', 'V7', 'iii7', 'vi7']),
        P('Diatonic seventh turnaround', 'Jazz', ['Imaj7', 'iii7', 'vi7', 'ii7', 'V7']),
        P('I–vi–ii–V7', 'Pop, jazz', ['I', 'vi', 'ii', 'V7']),
        P('Minor plagal "Amen"', 'Pop', ['I', 'iv', 'I']),
        P('Backdoor cadence in triads', 'Pop', ['I', 'iv', 'bVII', 'I']),
        P('Mixolydian vamp', 'Rock', ['I', 'bVII', 'I']),
        P('"Sweet Home Alabama"', 'Rock', ['I', 'bVII', 'IV']),
        P('Rock ♭VII–♭VI loop', 'Rock', ['I', 'bVII', 'bVI', 'bVII']),
        P('Minor rock cadence', 'Rock, film', ['i', 'bVI', 'bVII', 'i']),
        P('Minor loop through ♭III', 'Pop', ['i', 'bIII', 'bVII', 'iv']),
        P('Minor circle in triads', 'Pop', ['i', 'iv', 'bVII', 'bIII']),
        P('Minor with ♭VI', 'Classical, film', ['i', 'bVI', 'iv', 'V']),
        P('Minor tonic and subdominant', 'Classical', ['i', 'iv', 'i', 'V', 'i']),
        P('Baroque circle of fifths', 'Classical', ['I', 'IV', 'vii°', 'iii', 'vi', 'ii', 'V', 'I']),
      ],
    },
    hard: {
      name: 'Hard', desc: 'Secondary dominants, Neapolitan, jazz cycles', beat: 1.0,
      pad: ['I', 'i', 'bII', 'ii', 'ii°', 'II', 'bIII', 'iii', 'III', 'IV', 'iv', '#iv°', 'V', 'v', 'bVI', 'vi', 'VI', 'bVII', 'vii°'],
      list: [
        P('Turnaround with VI7', 'Jazz', ['Imaj7', 'VI7', 'ii7', 'V7']),
        P('Ragtime progression', 'Ragtime, jazz', ['I', 'VI7', 'II7', 'V7', 'I']),
        P('Circle to the tonic', 'Jazz', ['iii7', 'VI7', 'ii7', 'V7', 'Imaj7']),
        P('Dominant of the dominant', 'Classical', ['I', 'II7', 'V7', 'I']),
        P('Neapolitan cadence', 'Classical', ['i', 'bII', 'V7', 'i']),
        P('Secondary dominant to vi', 'Pop', ['I', 'III7', 'vi', 'IV']),
        P('Minor ii–V–i', 'Jazz', ['iiø7', 'V7', 'i']),
        P('Diminished passing chord', 'Classical', ['I', 'IV', '#iv°7', 'V']),
        P('Gospel plagal', 'Gospel, pop', ['I', 'I7', 'IV', 'iv', 'I']),
        P('Borrowed ♭III', 'Rock', ['I', 'bIII', 'IV', 'I']),
        P('Leading-tone cadence', 'Classical', ['I', 'ii', 'vii°7', 'I']),
        P('Minor circle of fifths', 'Jazz, baroque', ['i', 'iv7', 'bVII7', 'bIIImaj7', 'bVImaj7', 'iiø7', 'V7', 'i']),
        P('Major circle of fifths', 'Jazz, baroque', ['Imaj7', 'IVmaj7', 'iii7', 'vi7', 'ii7', 'V7', 'Imaj7']),
        P('Backdoor ii–V', 'Jazz', ['Imaj7', 'iv7', 'bVII7', 'Imaj7']),
        P('Double ii–V', 'Jazz', ['iii7', 'VI7', 'ii7', 'V7', 'I']),
        P('Rhythm changes, first bars', 'Jazz', ['Imaj7', 'VI7', 'ii7', 'V7', 'iii7', 'VI7', 'ii7', 'V7']),
        P('Tritone substitution', 'Jazz', ['ii7', 'bII7', 'Imaj7']),
        P('Chromatic descent with tritone subs', 'Jazz', ['iii7', 'bIII7', 'ii7', 'bII7', 'Imaj7']),
        P('Quick-change blues', 'Blues, jazz', ['I7', 'IV7', 'I7', 'I7', 'IV7', 'IV7', 'I7', 'I7', 'V7', 'IV7', 'I7', 'V7']),
        P('"Creep" progression', 'Rock', ['I', 'III', 'IV', 'iv']),
        P('Secondary dominant to IV', 'Pop', ['I', 'I7', 'IV', 'iv']),
        P('I–III7–IV–iv', 'Pop', ['I', 'III7', 'IV', 'iv']),
        P('Gospel turnaround', 'Gospel, blues', ['I', 'I7', 'IV', '#iv°7', 'I', 'V7', 'I']),
        P('Neapolitan after iv', 'Classical', ['i', 'iv', 'bII', 'V7', 'i']),
        P('Chromatic mediant', 'Film', ['I', 'bVI', 'I']),
        P('Picardy third', 'Classical, baroque', ['i', 'iv', 'V7', 'I']),
      ],
    },
  };
  for (const L of Object.values(LEVELS)) {
    for (const p of L.list) {
      p.parsed = p.chords.map(parse);
      p.minor = p.parsed[0].base === 'i' || p.parsed[p.parsed.length - 1].base === 'i';
      const missing = p.parsed.filter(c => !L.pad.includes(c.base));
      if (missing.length) console.warn(`${p.name}: ${missing.map(c => c.sym)} not on the ${L.name} pad`);
    }
  }
  const ALL_NUMERALS = LEVELS.hard.pad;
  const COUNTS = [5, 10, 20];
  const KEY = 'progressions.settings';

  const settings = App.loadSettings(KEY, { level: 'easy', count: 5, keyChord: true },
    s => s.level in LEVELS && COUNTS.includes(s.count));
  settings.keyChord = settings.keyChord !== false;
  settings.endless = settings.endless === true;
  const save = () => App.saveSettings(KEY, settings);

  // ---------- Voicing ----------
  // Tonics that keep every chord spelled with at most one sharp or flat
  const TONICS = [['C', 0], ['D', 0], ['E', -1], ['E', 0], ['F', 0], ['G', 0], ['A', 0], ['B', -1], ['A', -1], ['B', 0]];

  function spell(midi, L) {
    const alter = ((midi - SEMIS[L]) % 12 + 18) % 12 - 6;
    return { d: ((midi - alter - SEMIS[L]) / 12 - 1) * 7 + L, alter, midi };
  }

  // Bass plays the root; three upper voices move to the nearest notes of each next chord.
  function voice(tonic, parsed) {
    const [letter, alter] = tonic;
    const tl = LETTERS.indexOf(letter), tonicPc = (SEMIS[tl] + alter + 12) % 12;
    let prev = null;
    return parsed.map(c => {
      const rootL = (tl + c.degree) % 7;
      const rootPc = (tonicPc + MAJOR[c.degree] + c.shift + 12) % 12;
      const tones = c.semis.map((s, i) => ({ step: TONE_STEPS[i], pc: (rootPc + s) % 12 }));
      const upper = tones.length === 4 ? tones.slice(1) : tones;   // 7th chords: root stays in the bass
      let best = null;
      for (let first = 0; first < upper.length; first++) {
        for (let low = 52; low <= 70; low++) {
          if (low % 12 !== upper[first].pc) continue;
          const v = [{ midi: low, t: upper[first] }];
          for (let k = 1; k < upper.length; k++) {
            const t = upper[(first + k) % upper.length];
            let m = v[k - 1].midi + 1;
            while (m % 12 !== t.pc) m++;
            v.push({ midi: m, t });
          }
          if (v[v.length - 1].midi > 79) continue;
          const cost = prev ? v.reduce((a, x, i) => a + Math.abs(x.midi - prev[i]), 0) : Math.abs(v[1].midi - 67) * 2;
          if (!best || cost < best.cost) best = { cost, v };
        }
      }
      prev = best.v.map(x => x.midi);
      const bass = 40 + ((rootPc - 40) % 12 + 12) % 12;              // E2 to D♯3
      return [{ midi: bass, t: tones[0] }, ...best.v].map(({ midi, t }) => spell(midi, (rootL + t.step) % 7));
    });
  }

  function makeQuestion(L, prev) {
    for (;;) {
      const prog = pick(L.list);
      if (prev && prev.prog === prog && L.list.length > 1) continue;
      const tonic = pick(TONICS);
      const chords = voice(tonic, prog.parsed);
      if (chords.some(ch => ch.some(n => Math.abs(n.alter) > 1))) continue;
      const [letter, alter] = tonic;
      const keyName = `${letter}${alter === -1 ? '♭' : alter === 1 ? '♯' : ''} ${prog.minor ? 'minor' : 'major'}`;
      // Home chord for reference: I or i
      const home = voice(tonic, [parse(prog.minor ? 'i' : 'I')])[0];
      return { prog, chords, home, keyName };
    }
  }

  // ---------- Setup ----------
  const levelGrid = $('pLevelGrid'), countGroup = $('pCountGroup'), chordToggle = $('pKeyChordToggle');
  const syncCount = App.countPicker(countGroup, COUNTS, settings, () => { save(); syncSetup(); }, 'pcount');

  function buildSetup() {
    levelGrid.innerHTML = Object.entries(LEVELS).map(([id, L]) => `
      <button type="button" class="level-tile" role="radio" id="plevel-${id}" data-level="${id}" aria-checked="false">
        <span class="level-name">${L.name}</span>
        <span class="level-notes">${L.list.length} progressions</span>
        <span class="level-desc">${L.desc}</span>
      </button>`).join('');
    levelGrid.addEventListener('click', e => {
      const b = e.target.closest('[data-level]');
      if (!b) return;
      settings.level = b.dataset.level;
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
    syncCount();
    chordToggle.setAttribute('aria-checked', String(settings.keyChord));
  }

  // ---------- Game ----------
  const stage = $('pStage'), card = $('pCard'), scroll = $('pScroll'), svg = $('pStaff');
  const caption = $('pCaption'), pad = $('pPad'), confirmBtn = $('pConfirm');
  const timers = App.makeTimers();
  let game = null;

  function playChord(notes, at, dur) {
    const vol = 1.2 / Math.sqrt(notes.length);
    for (const n of notes) App.play(n.midi, at, vol, dur);
  }

  function playProgression(withHome) {
    const { q, level } = game;
    App.stopAll();
    let t = 0;
    if (withHome && settings.keyChord) {
      playChord(q.home, 0, 1.0);
      t = 1.7;
    }
    q.chords.forEach((notes, i) => playChord(notes, t + i * level.beat, level.beat * 0.92));
  }

  // The key's home chord on its own, as a reference
  function playHome() {
    App.ensureAudio();
    if (!game || !game.q) return;
    App.stopAll();
    playChord(game.q.home, 0, 1.4);
  }

  function buildPad(L) {
    // One column per scale degree, variants stacked (♭II, ii, ii°, II …)
    const cols = ROMANS.map(() => []);
    for (const id of L.pad) cols[parse(id).degree].push(id);
    pad.innerHTML = cols.map(ids => `<div class="pad-col">${ids.map(id =>
      `<button type="button" class="pad-btn" data-num="${id}">${pretty(id)}</button>`).join('')}</div>`).join('');
  }

  async function start() {
    const loading = App.ensureAudio();
    timers.clear();
    const level = LEVELS[settings.level];
    const endless = settings.endless;
    game = { level, endless, lives: App.LIVES, i: 0, total: endless ? Infinity : settings.count, correct: 0, chords: 0, perfect: 0, times: [], q: null, phase: 'loading', input: [] };
    buildPad(level);
    App.show('progGame');
    caption.textContent = 'Loading piano…';
    showProgress();
    update();
    const g = game;
    await Promise.race([loading, App.wait(5000)]);
    if (game !== g) return;   // left the screen while loading
    next();
  }

  function next() {
    const q = makeQuestion(game.level, game.q);
    Object.assign(game, { q, input: [], phase: 'input', start: 0 });
    caption.textContent = `${q.keyName.split(' ')[1] === 'minor' ? 'Minor' : 'Major'} key · ${q.chords.length} chords`;
    showProgress();
    update();
    draw(true);
    timers.later(() => {
      if (!game) return;
      playProgression(true);
      game.start = performance.now();
    }, 300);
  }

  function showProgress() {
    if (game.endless) App.setRun($('pProgressFill'), $('pProgressText'), game.lives, game.perfect);
    else App.setProgress($('pProgressFill'), $('pProgressText'), game.i, game.total);
  }
  const isLast = () => (game.endless ? game.lives <= 0 : game.i + 1 >= game.total);

  function update() {
    const n = game.q ? game.q.chords.length : 0;
    $('pHome').disabled = !game.q;
    $('pDelete').disabled = game.phase !== 'input' || !game.input.length;
    for (const b of pad.querySelectorAll('.pad-btn')) b.disabled = game.phase !== 'input';
    if (game.phase === 'review') {
      confirmBtn.textContent = isLast() ? 'See results' : 'Next progression';
      confirmBtn.disabled = false;
    } else {
      confirmBtn.textContent = 'Confirm';
      confirmBtn.disabled = game.phase !== 'input' || game.input.length < n;
    }
  }

  function draw(instant) {
    if (!game || !game.q) return;
    const { q, input, phase } = game;
    const n = q.chords.length;
    const cols = q.chords.map((notes, i) => {
      if (phase !== 'review') {
        return { notes: [], cursor: i === input.length, label: input[i] ? pretty(input[i]) : '' };
      }
      const truth = q.prog.parsed[i];
      const ok = input[i] === truth.base;
      const col = { notes: notes.map(x => ({ ...x, state: ok ? 'ok' : 'bad' })), label: pretty(input[i]), labelState: ok ? 'ok' : 'bad' };
      if (!ok) Object.assign(col, { sublabel: pretty(truth.sym), sublabelState: 'ok' });
      return col;
    });
    const layout = App.fitInStage(stage, card, scroll, maxHeight =>
      App.renderGrandStaff(svg, cols, { cols: n, colWidth: 78, sublabels: true, fit: scroll.clientWidth, min: 0.7, max: 1.4, maxHeight }));
    const focus = phase === 'review' ? 0 : Math.min(input.length, n - 1);
    scroll.scrollTo({
      left: Math.max(0, layout.centers[focus] * layout.scale - scroll.clientWidth * 0.55),
      behavior: instant || reduceMotion.matches ? 'auto' : 'smooth',
    });
  }

  function enter(id) {
    App.ensureAudio();
    if (!game || game.phase !== 'input' || game.input.length >= game.q.chords.length) return;
    game.input.push(id);
    update();
    draw();
  }

  function removeLast() {
    if (!game || game.phase !== 'input' || !game.input.length) return;
    game.input.pop();
    update();
    draw();
  }

  function confirm() {
    if (!game || !game.q) return;
    if (game.phase === 'review') return advance();
    const { q, input } = game;
    const n = q.chords.length;
    if (game.phase !== 'input' || input.length < n) return;
    game.phase = 'review';
    const ms = performance.now() - (game.start || performance.now());
    game.times.push(ms);
    let right = 0;
    q.prog.parsed.forEach((c, i) => {
      const ok = input[i] === c.base;
      if (ok) right++;
      App.recordItem('progressions', c.base, ok, ms / n, { label: pretty(c.base) });
      if (!ok) App.recordMixup('progressions', c.base, input[i]);
    });
    game.correct += right;
    game.chords += n;
    if (right === n) game.perfect++;
    else game.lives--;
    if (game.endless) showProgress();
    caption.textContent = `${q.prog.name} · ${q.prog.style} · ${q.keyName}`;
    update();
    draw();
    if (right < n) timers.later(() => game && playProgression(false), 400);   // hear it again
  }

  function advance() {
    const last = isLast();
    game.i++;
    if (last) return finish();
    next();
  }

  function finish() {
    const g = game;
    App.showResults({
      correct: g.endless ? g.perfect : g.correct,
      total: g.endless ? g.i : g.chords,
      times: g.times,
      unit: 'progression',
      statsKey: 'progressions',
      endless: g.endless ? { level: settings.level, name: g.level.name } : null,
      summary: g.endless
        ? `${g.perfect} perfect ${g.perfect === 1 ? 'progression' : 'progressions'} before ${App.LIVES} misses · ${g.correct} of ${g.chords} chords correct.`
        : `${g.correct} of ${g.chords} chords correct · ${g.perfect} of ${g.total} progressions perfect.`,
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
    App.show('progSetup');
  }

  // ---------- Wiring ----------
  buildSetup();
  syncSetup();
  $('pStartBtn').addEventListener('click', start);
  $('pQuitBtn').addEventListener('click', openSetup);
  $('pReplay').addEventListener('click', () => { App.ensureAudio(); if (game && game.q) playProgression(false); });
  $('pHome').addEventListener('click', playHome);
  $('pDelete').addEventListener('click', removeLast);
  confirmBtn.addEventListener('click', confirm);
  pad.addEventListener('click', e => {
    const b = e.target.closest('[data-num]');
    if (b) enter(b.dataset.num);
  });
  document.addEventListener('keydown', e => {
    if (e.repeat || e.metaKey || e.ctrlKey || e.altKey) return;
    if (App.current() !== 'progGame') return;
    if (e.key === 'Escape') return openSetup();
    if (e.key === 'Enter') { e.preventDefault(); return confirm(); }
    if (e.key === 'Backspace') { e.preventDefault(); return removeLast(); }
    if (e.key === ' ') { e.preventDefault(); if (game && game.q) playProgression(false); }
    if (e.key === 'h' || e.key === 'H') playHome();
  });

  App.watchStage(stage, () => { if (App.current() === 'progGame') draw(true); });
  App.onRedraw(() => { if (App.current() === 'progGame') draw(true); });

  window.Progressions = {
    open: openSetup,
    // For the stats page: every numeral, and its chord built in C major
    order: ALL_NUMERALS,
    label: pretty,
    shape: id => {
      // Root position, stacked from the root above middle C
      const c = parse(id), rootMidi = 60 + MAJOR[c.degree] + c.shift;
      return c.semis.map((s, i) => spell(rootMidi + s, (c.degree + TONE_STEPS[i]) % 7));
    },
  };
})();
