// Hub: pick a drill.
(() => {
  'use strict';
  const { $ } = App;

  $('openClefs').addEventListener('click', () => Clefs.open());
  $('openIntervals').addEventListener('click', () => Intervals.open());
  $('openChords').addEventListener('click', () => Chords.open());
  $('openMelody').addEventListener('click', () => Melody.open());
  $('openProgressions').addEventListener('click', () => Progressions.open());
  $('openIntonation').addEventListener('click', () => Intonation.open());
  $('openPitch').addEventListener('click', () => Pitch.open());
  $('openStandards').addEventListener('click', () => Standards.open());
  $('openStats').addEventListener('click', () => Stats.open());
  for (const b of document.querySelectorAll('[data-hub]')) b.addEventListener('click', () => App.show('hub'));

  App.onRedraw(() => {
    const box = $('hubClefsArt').parentElement.clientWidth;
    const opts = { barline: false, range: [11.6, -4], fit: box - 6, min: 0.2, max: 0.5 };
    App.renderStaff($('hubClefsArt'), [{ type: 'clef', clef: 'treble' }], opts);
    App.renderStaff($('hubIntervalsArt'), [
      { type: 'clef', clef: 'treble' },
      { type: 'note', p: 2, index: 0 },
      { type: 'note', p: 6, index: 1 },
    ], opts);
    App.renderStaff($('hubChordsArt'), [
      { type: 'clef', clef: 'treble' },
      { type: 'chord', ps: [2, 4, 6] },
      { type: 'chord', ps: [3, 5, 7] },
    ], opts);
    App.renderStaff($('hubProgArt'), [
      { type: 'clef', clef: 'treble' },
      { type: 'chord', ps: [0, 2, 4] },
      { type: 'chord', ps: [1, 3, 5] },
      { type: 'chord', ps: [2, 4, 6] },
    ], opts);
    App.renderStaff($('hubIntonationArt'), [
      { type: 'clef', clef: 'treble' },
      { type: 'note', p: 2, index: 0 },
      { type: 'note', p: 6, acc: 'sharp', index: 1 },
    ], opts);
    App.renderStaff($('hubMelodyArt'), [
      { type: 'clef', clef: 'treble' },
      { type: 'note', p: 2, index: 0 },
      { type: 'note', p: 5, index: 1 },
      { type: 'note', p: 4, index: 2 },
    ], opts);
  });

  App.boot();
})();
