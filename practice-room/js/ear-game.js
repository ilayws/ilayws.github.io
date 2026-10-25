// Shared game screen for the listening drills (Intervals, Chords):
// play a question, show its first note on a grand staff, pick an answer, review mistakes.
window.EarGame = (() => {
  'use strict';
  const { $ } = App;

  const RIGHT_PAUSE = 1300;              // after a wrong answer, the player presses Continue
  const stage = $('earStage'), card = $('earCard');
  const grandSvg = $('grand'), grandWrap = $('grandWrap');
  const answerGrid = $('answerGrid'), caption = $('earCaption');
  const continueBtn = $('continueBtn'), stageRow = $('earStageRow'), extraBtn = $('earExtraBtn');
  const timers = App.makeTimers();

  /**
   * round: {
   *   total, unit ('interval' | 'chord'), statsKey ('intervals' | 'chords'), colWidth,
   *   choices: [{ id, label, name }],
   *   makeQuestion(prev) -> { answerId, ... },
   *   play(q), caption(q),
   *   columns(q, guessId | null) -> 3 grand-staff columns (q.answer is set once answered),
   *   onAgain(), onSettings()
   *   optional: sublabels (room for a second label line), rightPause (ms),
   *   statsItem(q) -> { id, label } (what to track instead of the answer),
   *   noMixups, describeMisses(misses) -> text for the results screen,
   *   extra: { label, play(q) } (a second sound button once answered),
   *   keys: { key: answerId } (computer-keyboard shortcuts)
   * }
   */
  let round = null;
  let game = null;

  // Cut off whatever is still sounding so replays never overlap
  function playFresh(q) {
    App.stopAll();
    round.play(q);
  }

  function staffOpts() {
    return {
      cols: 3,
      colWidth: round.colWidth,
      fit: grandWrap.clientWidth - 16,
      min: 0.5,
      max: 1.7,
      sublabels: !!round.sublabels,
    };
  }

  function buildAnswers() {
    answerGrid.innerHTML = round.choices.map(c =>
      `<button type="button" class="answer-btn" id="ans-${c.id}" data-ans="${c.id}" aria-label="${c.name}" title="${c.name}">${c.label}</button>`).join('');
    answerGrid.dataset.unit = round.unit;
    answerGrid.dataset.kind = round.statsKey;
  }

  // Balanced rows: 8 buttons -> 4 + 4, 13 -> 7 + 6, 25 -> 7 + 7 + 7 + 4
  function layoutAnswers() {
    const n = answerGrid.children.length;
    if (!n) return;
    const narrow = answerGrid.clientWidth < 480;
    const maxCols = narrow ? (round.unit === 'chord' ? 6 : 5) : 7;
    const rows = Math.ceil(n / maxCols);
    answerGrid.style.setProperty('--cols', Math.ceil(n / rows));
  }

  function setContinue(visible) {
    continueBtn.hidden = !visible;
    stageRow.classList.toggle('has-continue', visible);
  }

  async function start(r) {
    round = r;
    const loading = App.ensureAudio();
    timers.clear();
    buildAnswers();
    setContinue(false);
    extraBtn.hidden = true;
    game = { q: null, i: 0, total: r.total, results: [], times: [], misses: [], locked: true };
    App.show('earGame');
    layoutAnswers();
    caption.textContent = 'Loading piano…';
    App.setProgress($('earProgressFill'), $('earProgressText'), 0, game.total);
    render([null, null, null]);
    const g = game;
    await Promise.race([loading, App.wait(5000)]);
    if (game !== g) return;   // left the screen while loading
    nextQuestion();
  }

  function nextQuestion() {
    setContinue(false);
    extraBtn.hidden = true;
    game.q = round.makeQuestion(game.q);
    game.q.answer = null;
    game.q.start = 0;
    game.locked = false;
    for (const b of answerGrid.children) b.classList.remove('ok', 'bad');
    caption.textContent = round.caption(game.q);
    App.setProgress($('earProgressFill'), $('earProgressText'), game.i, game.total);
    draw();
    timers.later(() => {
      if (!game) return;
      playFresh(game.q);
      game.q.start = performance.now();
    }, 250);
  }

  function render(columns) {
    App.fitInStage(stage, card, grandWrap, maxHeight =>
      App.renderGrandStaff(grandSvg, columns, { ...staffOpts(), maxHeight }));
  }

  function draw() {
    if (!game || !game.q) return;
    const a = game.q.answer;
    render(round.columns(game.q, a && !a.ok ? a.id : null));
  }

  function answer(id) {
    App.ensureAudio();
    if (App.current() !== 'earGame' || !game || game.locked || !game.q.start) return;
    const q = game.q;
    const ok = id === q.answerId;
    game.locked = true;
    const ms = performance.now() - q.start;
    game.times.push(ms);
    game.results.push(ok);
    const label = id => round.choices.find(c => c.id === id).label;
    const item = round.statsItem ? round.statsItem(q) : { id: q.answerId, label: label(q.answerId) };
    App.recordItem(round.statsKey, item.id, ok, ms, { label: item.label });
    if (!ok && !round.noMixups) App.recordMixup(round.statsKey, q.answerId, id);
    if (round.extra) { extraBtn.textContent = round.extra.label; extraBtn.hidden = false; }
    q.answer = { id, ok };
    $(`ans-${q.answerId}`).classList.add('ok');
    if (!ok) {
      $(`ans-${id}`).classList.add('bad');
      game.misses.push([q.answerId, id]);
      timers.later(() => game && playFresh(q), 300);   // hear the right one again
    }
    draw();
    if (ok) {
      timers.later(advance, round.rightPause || RIGHT_PAUSE);
    } else {
      setContinue(true);
      continueBtn.focus({ preventScroll: true });
    }
  }

  function continueAfterMistake() {
    if (!game || continueBtn.hidden) return;
    setContinue(false);
    advance();
  }

  function advance() {
    game.i++;
    if (game.i >= game.total) return finish();
    nextQuestion();
  }

  function finish() {
    const label = id => round.choices.find(c => c.id === id).label;
    const counts = new Map();
    for (const [right, said] of game.misses) {
      const k = `${label(right)} as ${label(said)}`;
      counts.set(k, (counts.get(k) || 0) + 1);
    }
    const mixed = [...counts].sort((a, b) => b[1] - a[1]).map(([k, n]) => (n > 1 ? `${k} (${n}×)` : k));
    App.showResults({
      correct: game.results.filter(Boolean).length,
      total: game.total,
      times: game.times,
      unit: round.unit,
      statsKey: round.statsKey,
      extra: round.describeMisses ? round.describeMisses(game.misses) : mixed.length ? `Mixed up: ${mixed.join(', ')}.` : '',
      review: null,
      onAgain: round.onAgain,
      onSettings: round.onSettings,
    });
    game = null;
  }

  function stop() {
    timers.clear();
    setContinue(false);
    extraBtn.hidden = true;
    game = null;
  }

  // ---------- Wiring ----------
  $('earQuitBtn').addEventListener('click', () => { stop(); round.onSettings(); });
  $('replayBtn').addEventListener('click', () => {
    App.ensureAudio();
    if (game && game.q && game.q.start) playFresh(game.q);
  });
  continueBtn.addEventListener('click', continueAfterMistake);
  extraBtn.addEventListener('click', () => {
    if (!game || !game.q || !round.extra) return;
    App.stopAll();
    round.extra.play(game.q);
  });
  answerGrid.addEventListener('click', e => {
    const b = e.target.closest('[data-ans]');
    if (b) answer(b.dataset.ans);
  });

  document.addEventListener('keydown', e => {
    if (e.repeat || e.metaKey || e.ctrlKey || e.altKey) return;
    if (App.current() !== 'earGame') return;
    if (e.key === 'Escape') { stop(); return round.onSettings(); }
    if (e.key === 'Enter' && !continueBtn.hidden) {
      e.preventDefault();
      return continueAfterMistake();
    }
    if (round && round.keys && round.keys[e.key] && game && !game.locked) {
      e.preventDefault();
      return answer(round.keys[e.key]);
    }
    if (e.key === ' ' && game && game.q && game.q.start) {
      e.preventDefault();
      playFresh(game.q);
    }
  });

  App.watchStage(stage, () => { if (App.current() === 'earGame' && round) draw(); });
  App.onRedraw(() => {
    if (App.current() !== 'earGame' || !round) return;
    layoutAnswers();
    draw();
  });

  return { start, stop };
})();
