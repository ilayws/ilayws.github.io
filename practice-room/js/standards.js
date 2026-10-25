// Standard Guesser: hear a jazz recording (audio only, from YouTube) and type the tune's name.
// Endless: three hearts, a wrong answer or a skip costs one. Faster answers score more.
(() => {
  'use strict';
  const { $ } = App;

  const SONGS = window.STANDARDS;
  const MAX_POINTS = 1000, MIN_POINTS = 100, HALF_LIFE = 15;   // a song's points halve every 15 s of listening
  const RIGHT_PAUSE = 2200;
  const STALL_MS = 5000;      // still silent this long after loading: the browser blocked autoplay
  const MAX_ERRORS = 4;       // songs in a row that won't play before giving up

  const points = secs => Math.max(MIN_POINTS, Math.round((MAX_POINTS * 0.5 ** (secs / HALF_LIFE)) / 10) * 10);
  const clock = s => `${Math.floor(s / 60)}:${String(Math.floor(s % 60)).padStart(2, '0')}`;
  const fmt = n => n.toLocaleString('en-US');

  // ---------- Lenient matching ----------
  // Compare without case, accents, punctuation, spaces or a leading article; allow about
  // one typo per four letters, extra words around the title, or most of a long title.
  // A guess that is closer to another tune in the list ("Blue Moon" for Blue Monk) doesn't count.
  const STOP = new Set(['the', 'a', 'an']);
  const words = s => s.normalize('NFKD').replace(/[̀-ͯ]/g, '').toLowerCase()
    .replace(/&/g, ' and ').replace(/['’`]/g, '').replace(/[^a-z0-9]+/g, ' ').trim().split(' ').filter(Boolean);

  function spellings(name) {
    const out = new Set();
    for (const n of [name, name.replace(/\(.*?\)/g, ' ')]) {
      const w = words(n);
      out.add(w.join(''));
      out.add(w.filter(x => !STOP.has(x)).join(''));
    }
    out.delete('');
    return [...out];
  }

  function distance(a, b) {
    let row = Array.from({ length: b.length + 1 }, (_, j) => j);
    for (let i = 1; i <= a.length; i++) {
      const next = [i];
      for (let j = 1; j <= b.length; j++) {
        next[j] = Math.min(row[j] + 1, next[j - 1] + 1, row[j - 1] + (a[i - 1] === b[j - 1] ? 0 : 1));
      }
      row = next;
    }
    return row[b.length];
  }

  for (const s of SONGS) s.keys = [s.t, ...(s.a || [])].flatMap(spellings);

  // How far a guess is from a song: typos, or 0.5 for the title inside the guess and 1 for
  // most of a long title; Infinity when it's too far to count.
  function cost(gs, song) {
    let best = Infinity;
    for (const k of song.keys) {
      // Names with digits ("Take 5") have to be exact: a typo there is a different name
      const allowed = k.length < 4 || /\d/.test(k) ? 0 : k.length < 8 ? 1 : Math.min(5, Math.floor(k.length / 4));
      for (const g of gs) {
        const d = distance(g, k);
        if (d <= allowed) best = Math.min(best, d);
        if (k.length >= 5 && g.includes(k)) best = Math.min(best, 0.5);                         // "autumn leaves by cannonball"
        if (g.length >= 6 && k.includes(g) && g.length >= 0.7 * k.length) best = Math.min(best, 1);   // most of a long title
      }
    }
    return best;
  }

  function matches(guess, song) {
    const gs = spellings(guess);
    const c = cost(gs, song);
    return c < Infinity && !SONGS.some(other => other !== song && cost(gs, other) < c);
  }

  // ---------- YouTube player (kept out of sight) ----------
  let apiReady = null, player = null, playerReady = null;

  function loadApi() {
    if (!apiReady) {
      apiReady = new Promise((resolve, reject) => {
        if (window.YT && window.YT.Player) return resolve();
        window.onYouTubeIframeAPIReady = resolve;
        const s = document.createElement('script');
        s.src = 'https://www.youtube.com/iframe_api';
        s.onerror = () => { apiReady = null; reject(new Error('offline')); };
        document.head.appendChild(s);
      });
    }
    return apiReady;
  }

  function makePlayer() {
    if (!playerReady) {
      playerReady = new Promise(resolve => {
        player = new YT.Player('sgPlayer', {
          width: 200, height: 112,
          playerVars: { autoplay: 1, controls: 0, disablekb: 1, fs: 0, iv_load_policy: 3, playsinline: 1, rel: 0 },
          events: { onReady: () => resolve(), onStateChange: onState, onError: onPlayError },
        });
      });
    }
    return playerReady;
  }

  // ---------- Game ----------
  const card = $('sgCard'), input = $('sgInput'), guessBtn = $('sgGuess'), skipBtn = $('sgSkip'), nextBtn = $('sgNext');
  const timers = App.makeTimers();
  let game = null, ticker = 0;

  // Songs come in a shuffled order and don't repeat until every one has played
  function draw() {
    if (!game.bag.length) {
      game.bag = SONGS.filter(s => !game.broken.has(s) && s !== game.song);
      for (let i = game.bag.length - 1; i > 0; i--) {
        const j = Math.floor(Math.random() * (i + 1));
        [game.bag[i], game.bag[j]] = [game.bag[j], game.bag[i]];
      }
    }
    return game.bag.pop();
  }

  // Seconds of the song heard so far (the clock stops while it buffers)
  function elapsed() {
    if (game.t0 == null) return 0;
    return game.looped + Math.max(0, player.getCurrentTime() - game.t0);
  }

  function showProgress() {
    App.setRun($('sgProgressFill'), $('sgProgressText'), game.lives, game.songs, `${fmt(game.score)} pts`);
  }

  function setPanel(name) {
    card.dataset.panel = name;   // 'now' | 'tap' | 'reveal' | 'error'
  }

  function setStatus(text) {
    $('sgStatus').textContent = text;
  }

  function updateClock() {
    const t = elapsed();
    $('sgClock').textContent = clock(t);
    $('sgWorth').textContent = `Worth ${fmt(points(t))}`;
  }

  function setControls() {
    const listening = game.phase === 'listening';
    guessBtn.disabled = !listening;
    skipBtn.disabled = !listening;
    input.disabled = game.phase === 'reveal';
    nextBtn.hidden = game.phase !== 'reveal';
    skipBtn.hidden = game.phase === 'reveal';
    nextBtn.textContent = game.lives <= 0 ? 'See results' : 'Next song';
  }

  async function start() {
    timers.clear();
    const g = game = {
      lives: App.LIVES, score: 0, songs: 0, guessMs: [], missed: [],
      bag: [], broken: new Set(), song: null, phase: 'loading', t0: null, looped: 0, errors: 0,
    };
    App.show('standardsGame');
    showProgress();
    setPanel('now');
    setStatus('Loading the player…');
    $('sgClock').textContent = '0:00';
    $('sgWorth').textContent = `Worth ${fmt(MAX_POINTS)}`;
    input.value = '';
    setControls();
    try {
      await loadApi();
      await makePlayer();
    } catch (e) {
      if (game === g) showError();
      return;
    }
    if (game !== g) return;   // left while the player loaded
    nextSong();
  }

  function nextSong() {
    timers.clear();
    Object.assign(game, { song: draw(), phase: 'loading', t0: null, looped: 0 });
    setPanel('now');
    setStatus('Loading…');
    card.classList.remove('playing');
    updateClock();
    input.value = '';
    setControls();
    player.loadVideoById({ videoId: game.song.v, startSeconds: 0 });
    timers.later(() => {
      if (game && game.phase === 'loading') {
        setPanel('tap');   // the browser wants a tap on the player itself before it plays sound
        setStatus('Tap ▶ to start the music');
      }
    }, STALL_MS);
  }

  function onState(e) {
    if (!game || App.current() !== 'standardsGame') return;
    const S = YT.PlayerState;
    card.classList.toggle('playing', e.data === S.PLAYING);
    if (e.data === S.PLAYING) {
      game.errors = 0;
      if (game.phase === 'loading') {
        game.phase = 'listening';
        game.t0 = player.getCurrentTime();
        setPanel('now');
        setStatus('Name this tune');
        setControls();
        input.focus({ preventScroll: true });
        clearInterval(ticker);
        ticker = setInterval(() => { if (game && game.phase === 'listening') updateClock(); }, 250);
      }
    } else if (e.data === S.ENDED && game.phase !== 'over') {
      // Keep it going until there's an answer
      if (game.t0 != null) { game.looped += Math.max(0, player.getDuration() - game.t0); game.t0 = 0; }
      player.seekTo(0, true);
      player.playVideo();
    }
  }

  // A video that won't play here: quietly move on, without costing a heart
  function onPlayError() {
    if (!game || App.current() !== 'standardsGame' || game.phase === 'reveal' || game.phase === 'over') return;
    game.broken.add(game.song);
    if (++game.errors >= MAX_ERRORS) return showError();
    nextSong();
  }

  function showError() {
    game.phase = 'over';
    timers.clear();
    setControls();
    guessBtn.disabled = skipBtn.disabled = input.disabled = true;
    $('sgErrorText').textContent = location.protocol === 'file:'
      ? 'YouTube recordings can’t play from a file saved on this device. Open the online Practice Room instead:'
      : 'The recordings couldn’t be loaded from YouTube. Check your internet connection and try again.';
    $('sgErrorLink').hidden = location.protocol !== 'file:';
    setPanel('error');
  }

  function guess(e) {
    e.preventDefault();
    if (!game) return;
    if (game.phase === 'reveal') return next();
    const text = input.value.trim();
    if (game.phase !== 'listening' || !text) return;
    reveal(matches(text, game.song) ? 'right' : 'wrong', text);
  }

  function skip() {
    if (game && game.phase === 'listening') reveal('skip', '');
  }

  function reveal(kind, said) {
    clearInterval(ticker);
    const t = elapsed();
    game.phase = 'reveal';
    const s = game.song;
    let verdict;
    if (kind === 'right') {
      const pts = points(t);
      game.score += pts;
      game.songs++;
      game.guessMs.push(t * 1000);
      verdict = `Correct · +${fmt(pts)}`;
    } else {
      game.lives--;
      game.missed.push(s.t);
      verdict = kind === 'skip' ? 'Skipped' : 'Not quite';
    }
    $('sgVerdict').textContent = verdict;
    $('sgVerdict').className = `sg-verdict ${kind === 'right' ? 'ok' : 'bad'}`;
    $('sgTitle').textContent = s.t;
    $('sgBy').textContent = `${s.by} · named after ${clock(t)}`;
    $('sgSaid').textContent = kind === 'wrong' ? `You answered “${said}”` : '';
    $('sgSaid').hidden = kind !== 'wrong';
    $('sgLink').href = `https://www.youtube.com/watch?v=${s.v}`;
    setPanel('reveal');
    showProgress();
    setControls();
    if (kind === 'right') timers.later(next, RIGHT_PAUSE);
    else nextBtn.focus({ preventScroll: true });
  }

  function next() {
    if (!game || game.phase !== 'reveal') return;
    if (game.lives <= 0) return finish();
    nextSong();
  }

  function finish() {
    const g = game;
    g.phase = 'over';
    stopPlayer();
    const avgMs = g.guessMs.length ? g.guessMs.reduce((a, b) => a + b, 0) / g.guessMs.length : 0;
    const prev = App.recordGame('standards', { score: g.score, songs: g.songs, avgMs: Math.round(avgMs) });
    const records = [];
    if (g.score > prev.score) records.push('New high score!');
    if (g.songs > prev.songs) records.push(g.score > prev.score ? 'Most songs yet, too.' : 'Most songs yet!');
    const most = Math.max(prev.songs, g.songs);
    const best = `Best: ${fmt(Math.max(prev.score, g.score))} points, ${most} ${most === 1 ? 'song' : 'songs'}.`;
    const speed = g.songs ? ` ${(avgMs / 1000).toFixed(1)} s on average to name a song.` : '';
    App.showResults({
      custom: {
        title: 'Game over',
        tiles: [[fmt(g.score), 'Points'], [String(g.songs), g.songs === 1 ? 'Song named' : 'Songs named']],
        summary: `${records.join(' ')}${records.length ? ' ' : ''}${best}${speed}`,
      },
      correct: g.songs,
      total: g.songs + g.missed.length,
      times: g.guessMs,
      unit: 'song',
      extra: g.missed.length ? `Missed: ${g.missed.join(', ')}.` : '',
      review: null,
      onAgain: start,
      onSettings: openSetup,
    });
    game = null;
  }

  function stopPlayer() {
    timers.clear();
    clearInterval(ticker);
    card.classList.remove('playing');
    if (player && player.pauseVideo) player.pauseVideo();
  }

  function openSetup() {
    stopPlayer();
    game = null;
    const games = App.readStats().standards.games;
    $('sgBestScore').textContent = games.length ? fmt(Math.max(...games.map(x => x.score))) : '–';
    $('sgBestSongs').textContent = games.length ? String(Math.max(...games.map(x => x.songs))) : '–';
    App.show('standardsSetup');
  }

  // ---------- Wiring ----------
  $('sgCount').textContent = String(SONGS.length);
  $('sgStartBtn').addEventListener('click', start);
  $('sgQuitBtn').addEventListener('click', openSetup);
  $('sgForm').addEventListener('submit', guess);
  skipBtn.addEventListener('click', skip);
  nextBtn.addEventListener('click', next);
  $('sgErrorBack').addEventListener('click', openSetup);
  document.addEventListener('keydown', e => {
    if (App.current() !== 'standardsGame' || e.metaKey || e.ctrlKey || e.altKey) return;
    if (e.key === 'Escape') return openSetup();
    if (e.key === 'Enter' && game && game.phase === 'reveal' && document.activeElement !== input) { e.preventDefault(); next(); }
  });
  // Leaving the game screen any other way (back to the hub, the results) stops the music
  App.onShow((name, prev) => { if (prev === 'standardsGame' && name !== 'standardsGame') { stopPlayer(); if (game && game.phase !== 'over') game = null; } });

  window.Standards = { open: openSetup, matches };
})();
