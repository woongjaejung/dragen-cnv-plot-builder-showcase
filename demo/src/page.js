
    // 24 curated hues, then golden-angle HSL so any sample count stays distinguishable.
    const BASE_COLORS = ['#1f5f8b','#d17a22','#2c8a63','#b14d4d','#6b5b95','#118a99','#8e6c2f','#4b648d',
                         '#c43d8c','#3d9ac4','#8ab94a','#a14d8b','#5f7a3a','#bb6677','#0f766e','#a3421f',
                         '#2563a0','#7a7a1f','#a8326e','#356859','#9c5c1b','#4c4cae','#127a5c','#96331f'];
    const NTC_COLOR = '#8a8275';
    function paletteFor(n) {
      const out = BASE_COLORS.slice(0, Math.min(n, BASE_COLORS.length));
      for (let i = BASE_COLORS.length; i < n; i++) {
        const hue = (i * 137.508) % 360;
        const light = 32 + ((i % 3) * 9);
        out.push(`hsl(${hue.toFixed(1)},58%,${light}%)`);
      }
      return out;
    }
    const SYMBOLS = ['circle','square','diamond','triangle-up','cross','triangle-down','star','hexagon'];

    // scattergl is much faster for big genes at high sample counts, but some
    // managed Windows/VDI images ship without WebGL. Detect once and fall back
    // to SVG scatter rather than rendering an error message into the plot.
    const HAS_WEBGL = (() => {
      try {
        const c = document.createElement('canvas');
        return !!(c.getContext('webgl') || c.getContext('experimental-webgl'));
      } catch (e) { return false; }
    })();
    const GL_POINT_THRESHOLD = 6000;

    const FILTER_DESCRIPTIONS = {
      PASS: 'Passes all enabled CNV filters.',
      cnvQual: 'CNV QUAL is below the workflow threshold.',
      cnvLength: 'Segment is shorter than the minimum reportable target span (single noisy target).',
      cnvLikelihoodRatio: 'Log10 likelihood ratio of ALT to REF is below the workflow threshold.',
      dinucQual: 'Segment dinucleotide composition falls outside typical ranges, so the call is likely a false positive.',
    };
    const REFERENCE_LINES = {
      normalized_ratio: [
        { y: 1.0, color: 'rgba(46,120,72,0.55)', dash: 'solid', label: 'diploid 1.0' },
        { y: 0.5, color: 'rgba(177,77,77,0.55)', dash: 'dash',  label: 'het del 0.5' },
        { y: 1.5, color: 'rgba(31,95,139,0.55)', dash: 'dash',  label: 'het dup 1.5' },
      ],
      tn_log2: [
        { y: 0.0,   color: 'rgba(46,120,72,0.55)', dash: 'solid', label: 'diploid log2=0' },
        { y: -1.0,  color: 'rgba(177,77,77,0.55)', dash: 'dash',  label: 'het del log2=-1' },
        { y: 0.585, color: 'rgba(31,95,139,0.55)', dash: 'dash',  label: 'het dup log2=0.585' },
      ],
      segment_mean: [
        { y: 0.0, color: 'rgba(46,120,72,0.55)', dash: 'solid', label: 'segment 0' },
      ],
    };
    // Kept light on purpose: the band marks a region, the points carry the data.
    // Regions are merged before drawing (see buildShapes), so these alphas are
    // what actually lands on screen regardless of how many samples share a call.
    const BAND_COLORS = {
      DEL:  'rgba(177,77,77,0.11)',
      DUP:  'rgba(31,95,139,0.11)',
      LOSS: 'rgba(177,77,77,0.11)',
      GAIN: 'rgba(31,95,139,0.11)',
    };
    const BAND_COLOR_OTHER = 'rgba(110,104,95,0.10)';

    const meta = JSON.parse(document.getElementById('plot-meta').textContent);
    const ntcSet = new Set(meta.ntc_samples || []);
    const groupIndex = new Map(meta.groups.map((g, i) => [g, i]));
    const groupCache = new Map();

    // Each group lives in its own <script> tag and is parsed on first view, so
    // no single string ever approaches the JS max string length.
    function groupData(name) {
      const idx = groupIndex.get(name);
      if (idx === undefined) return null;
      if (groupCache.has(idx)) return groupCache.get(idx);
      const node = document.getElementById('g' + idx);
      if (!node) return null;
      const parsed = JSON.parse(node.textContent);
      groupCache.set(idx, parsed);
      return parsed;
    }

    function isDarkTheme() {
      return document.documentElement.getAttribute('data-theme') === 'dark';
    }
    function themePlotColors() {
      if (isDarkTheme()) {
        return {
          plotBg: '#0b1629',
          fontColor: '#e2ecf8',
          mutedColor: '#8ea6c7',
          gridColor: 'rgba(133,157,185,0.16)',
          zeroLineColor: 'rgba(133,157,185,0.34)',
          exonSepColor: 'rgba(133,157,185,0.28)',
          hoverBg: '#071a38',
          hoverBorder: '#00b5e2',
          hoverFont: '#edf4fc',
        };
      }
      return {
        plotBg: '#f8fafd',
        fontColor: '#0f1f38',
        mutedColor: '#61748f',
        gridColor: 'rgba(97,116,143,0.16)',
        zeroLineColor: 'rgba(97,116,143,0.35)',
        exonSepColor: 'rgba(97,116,143,0.30)',
        hoverBg: '#071a38',
        hoverBorder: '#00a3e0',
        hoverFont: '#ffffff',
      };
    }

    const themeLightBtn = document.getElementById('theme-light-btn');
    const themeDarkBtn = document.getElementById('theme-dark-btn');
    function applyTheme(theme, shouldRerender) {
      const resolved = theme === 'dark' ? 'dark' : 'light';
      document.documentElement.setAttribute('data-theme', resolved);
      if (themeLightBtn) themeLightBtn.classList.toggle('is-active', resolved === 'light');
      if (themeDarkBtn) themeDarkBtn.classList.toggle('is-active', resolved === 'dark');
      try { localStorage.setItem('dragen-cnv-report-theme', resolved); } catch (e) {}
      if (shouldRerender && window.__cnvReady) render(currentGroup());
    }
    if (themeLightBtn) themeLightBtn.addEventListener('click', () => applyTheme('light', true));
    if (themeDarkBtn) themeDarkBtn.addEventListener('click', () => applyTheme('dark', true));
    try {
      const savedTheme = localStorage.getItem('dragen-cnv-report-theme');
      if (savedTheme === 'dark' || savedTheme === 'light') applyTheme(savedTheme, false);
    } catch (e) {}

    const groupInput = document.getElementById('group-input');
    const groupList = document.getElementById('group-list');
    const metricSelect = document.getElementById('metric-select');
    const sortSelect = document.getElementById('sort-select');
    const renderBtn = document.getElementById('render-btn');
    const resetBtn = document.getElementById('reset-btn');
    const toggleRefs = document.getElementById('toggle-refs');
    const toggleBands = document.getElementById('toggle-bands');
    const toggleStrand = document.getElementById('toggle-strand');
    const toggleNtc = document.getElementById('toggle-ntc');
    const groupLabel = document.getElementById('group-label');
    const callsDesc = document.getElementById('calls-desc');
    const callsContent = document.getElementById('calls-content');
    const callsFilterNotes = document.getElementById('calls-filter-notes');
    const strandChip = document.getElementById('strand-chip');
    const ntcChip = document.getElementById('ntc-chip');

    // NTC samples all share one muted color, and do not consume a palette slot.
    const nonNtcCount = meta.samples.filter(s => !ntcSet.has(s)).length;
    const COLORS = paletteFor(nonNtcCount);
    const sampleColorByName = new Map();
    const sampleSymbolByIndex = new Map();
    let colorCursor = 0;
    meta.samples.forEach((sample, idx) => {
      if (ntcSet.has(sample)) {
        sampleColorByName.set(sample, NTC_COLOR);
        // Only the NTC styling itself supplies the X marker, so turning that
        // toggle off returns these samples to a plain circle.
        sampleSymbolByIndex.set(idx, 'circle');
      } else {
        sampleColorByName.set(sample, COLORS[colorCursor % COLORS.length]);
        sampleSymbolByIndex.set(idx, SYMBOLS[Math.floor(colorCursor / BASE_COLORS.length) % SYMBOLS.length]);
        colorCursor += 1;
      }
    });

    groupLabel.textContent = meta.group_label;
    document.getElementById('group-count-chip').textContent = `${meta.group_label} groups: ${meta.groups.length}`;
    document.getElementById('sample-count-chip').textContent = `samples: ${meta.samples.length}`;
    if (ntcSet.size) {
      ntcChip.hidden = false;
      ntcChip.textContent = `NTC: ${[...ntcSet].join(', ')}`;
    }
    document.getElementById('mode-chip').textContent = meta.has_gene_data ? 'mode: gene plot' : 'mode: contig plot';
    document.getElementById('dropped-chip').textContent = `unannotated rows skipped: ${meta.dropped_unannotated_rows || 0}`;

    // A datalist with tens of thousands of options is slow to build; cap it and
    // let the free-text box handle anything outside the shortlist.
    const DATALIST_CAP = 3000;
    for (const group of meta.groups.slice(0, DATALIST_CAP)) {
      const option = document.createElement('option');
      option.value = group;
      groupList.appendChild(option);
    }

    function metricLabel(metric) {
      if (metric === 'tn_log2') return 'Tangent-Normalized Log2';
      if (metric === 'segment_mean') return 'Segment Mean';
      return 'Normalized Ratio';
    }

    function sampleColor(sampleOrIndex) {
      const name = typeof sampleOrIndex === 'number' ? meta.samples[sampleOrIndex] : sampleOrIndex;
      return sampleColorByName.get(name) || BASE_COLORS[0];
    }

    // normalized_ratio is exactly 2**tn_log2, so it is derived here instead of shipped.
    function metricValue(data, i, metric) {
      if (metric === 'tn_log2') return data.v[i];
      if (metric === 'segment_mean') return data.m[i];
      const log2 = data.v[i];
      return log2 === null ? null : Math.pow(2, log2);
    }

    function filterBadgeClass(f) { return f === 'PASS' ? 'badge badge-pass' : 'badge badge-filtered'; }

    function renderFilterNotes(calls) {
      const filters = [...new Set(calls.map(c => c.call_filter || (c.status === 'PASS' ? 'PASS' : '')).filter(Boolean))];
      if (!filters.length) { callsFilterNotes.innerHTML = ''; return; }
      filters.sort((a, b) => (a === 'PASS' ? -1 : b === 'PASS' ? 1 : a.localeCompare(b)));
      const items = filters.map(f => {
        const d = FILTER_DESCRIPTIONS[f] || 'Meaning not annotated in this HTML.';
        return `<div class="filter-note"><span class="${filterBadgeClass(f)}">${f}</span><span>${d}</span></div>`;
      }).join('');
      callsFilterNotes.innerHTML = `<p class="filter-notes-title">Filter Definitions</p><div class="filter-notes-list">${items}</div>`;
    }

    function summarizeLabels(labels) {
      if (!labels || !labels.length) return '';
      if (!labels.every(l => /^-?\d+$/.test(String(l)))) return labels.join(', ');
      const nums = [...new Set(labels.map(l => parseInt(l, 10)))].sort((a, b) => a - b);
      const ranges = [];
      let start = nums[0], prev = nums[0];
      for (const n of nums.slice(1)) {
        if (n === prev + 1) { prev = n; continue; }
        ranges.push(start === prev ? `${start}` : `${start}-${prev}`);
        start = prev = n;
      }
      ranges.push(start === prev ? `${start}` : `${start}-${prev}`);
      return ranges.join(', ');
    }

    function hoverText(data, i, group, metric) {
      const rank = data.x[i];
      const t = data.targets[rank - 1];
      const value = metricValue(data, i, metric);
      const call = data.c[i] >= 0 ? data.calls[data.c[i]] : null;
      const sampleName = meta.samples[data.s[i]];
      const ntcTag = ntcSet.has(sampleName) ? ' [NTC]' : '';
      return [
        `sample: ${sampleName}${ntcTag}`,
        `${meta.group_label}: ${meta.has_gene_data ? group : t[1]}`,
        meta.has_gene_data ? `exon: ${t[0]}` : `x: ${t[0]}`,
        `plot target order: ${rank}`,
        `${metricLabel(metric)}: ${value === null || value === undefined ? 'NA' : value}`,
        `target: ${t[4]}`,
        `region: ${t[1]}:${t[2]}-${t[3]}`,
        call ? `call: ${call.call_type} CN=${call.call_copy_number || '?'} filter=${call.call_filter || ''} region=${call.contig}:${call.call_start}-${call.call_end}` : ''
      ].filter(Boolean).join('<br>');
    }

    function bandColorForCall(callType) {
      const key = String(callType || '').toUpperCase();
      if (BAND_COLORS[key]) return BAND_COLORS[key];
      if (key.includes('DEL') || key.includes('LOSS')) return BAND_COLORS.DEL;
      if (key.includes('DUP') || key.includes('GAIN')) return BAND_COLORS.DUP;
      return BAND_COLOR_OTHER;
    }

    function renderCalls(group, data) {
      const calls = data.calls;
      const axisHint = meta.has_gene_data
        ? 'Exons In Plot matches the x-axis exon labels above.'
        : 'Plot Target Range matches the target-order x-axis above.';
      callsDesc.textContent = `Rows below include PASS and non-PASS CNV calls whose regions overlap at least one displayed target in the selected ${meta.group_label}: ${group}. ${axisHint}`;
      if (!calls.length) {
        callsContent.innerHTML = '<div class="hint">No CNV calls overlapped this group.</div>';
        callsFilterNotes.innerHTML = '';
        return;
      }
      const exonHeader = meta.has_gene_data ? 'Exons In Plot' : 'Target Labels In Plot';
      let html = `<div class="calls-table-wrap"><table><thead><tr><th>Sample</th><th>Type</th><th>CN</th><th>Filter</th><th>Call Region</th><th>${exonHeader}</th><th>Plot Target Range</th></tr></thead><tbody>`;
      for (const call of calls) {
        const filterText = call.call_filter || (call.status === 'PASS' ? 'PASS' : '');
        const ntcBadge = ntcSet.has(call.sample_id) ? ' <span class="badge badge-ntc">NTC</span>' : '';
        const dot = `<span class="sample-chip"><span class="sample-dot" style="background:${sampleColor(call.sample_id)}"></span>${call.sample_id}${ntcBadge}</span>`;
        html += `<tr><td>${dot}</td><td>${call.call_type}</td><td>${call.call_copy_number || ''}</td><td><span class="${filterBadgeClass(filterText)}">${filterText}</span></td><td class="mono">${call.contig}:${call.call_start}-${call.call_end}</td><td>${summarizeLabels(call.exon_labels) || ''}</td><td class="mono">targets ${call.x_start ?? ''}-${call.x_end ?? ''} (${call.target_count})</td></tr>`;
      }
      callsContent.innerHTML = html + '</tbody></table></div>';
      renderFilterNotes(calls);
    }

    function buildShapes(data, metric) {
      const shapes = [];

      if (toggleBands.checked) {
        // One band per region, not per call. Calls are per sample, so a region
        // that is PASS in N samples used to stack N translucent rectangles; at
        // 30 samples a 0.16 alpha compounds to ~99% opaque and buries the
        // points it is supposed to sit behind. Merge the spans per color first.
        const spansByColor = new Map();
        for (const call of data.calls) {
          if ((call.call_filter || '').toUpperCase() !== 'PASS') continue;
          if (call.x_start == null || call.x_end == null) continue;
          const color = bandColorForCall(call.call_type);
          let spans = spansByColor.get(color);
          if (!spans) { spans = []; spansByColor.set(color, spans); }
          spans.push([call.x_start - 0.5, call.x_end + 0.5]);
        }
        const bandRect = (x0, x1, color) => ({
          type: 'rect', xref: 'x', yref: 'paper',
          x0, x1, y0: 0, y1: 1,
          line: { width: 0 }, fillcolor: color, layer: 'below',
        });
        for (const [color, spans] of spansByColor) {
          spans.sort((a, b) => a[0] - b[0]);
          let [from, to] = spans[0];
          for (let i = 1; i < spans.length; i++) {
            const [ns, ne] = spans[i];
            if (ns <= to) { if (ne > to) to = ne; continue; }
            shapes.push(bandRect(from, to, color));
            from = ns; to = ne;
          }
          shapes.push(bandRect(from, to, color));
        }
      }

      const tc = themePlotColors();
      for (const pos of (data.bp || [])) {
        shapes.push({
          type: 'line', xref: 'x', yref: 'paper',
          x0: pos, x1: pos, y0: 0, y1: 1,
          line: { color: tc.exonSepColor, width: 1, dash: 'dot' },
        });
      }

      if (toggleRefs.checked) {
        for (const ref of (REFERENCE_LINES[metric] || [])) {
          shapes.push({
            type: 'line', xref: 'paper', yref: 'y',
            x0: 0, x1: 1, y0: ref.y, y1: ref.y,
            line: { color: ref.color, width: 1.2, dash: ref.dash }, layer: 'below',
          });
        }
      }

      return shapes;
    }

    function buildAnnotations(metric) {
      if (!toggleRefs.checked) return [];
      return (REFERENCE_LINES[metric] || []).map(ref => ({
        xref: 'paper', yref: 'y',
        x: 1.005, y: ref.y, xanchor: 'left', yanchor: 'middle',
        text: ref.label, showarrow: false,
        font: { size: 10, color: ref.color },
      }));
    }

    function render(group) {
      const data = groupData(group);
      if (!data) { alert(`No data found for "${group}"`); return; }

      const metric = metricSelect.value;
      const sortMode = sortSelect.value;
      const bySample = new Map();

      for (let i = 0; i < data.s.length; i++) {
        const y = metricValue(data, i, metric);
        if (y === null || y === undefined) continue;
        const si = data.s[i];
        let trace = bySample.get(si);
        if (!trace) {
          const sampleName = meta.samples[si];
          const isNtc = ntcSet.has(sampleName);
          trace = {
            x: [], y: [], text: [], si, isNtc,
            name: isNtc ? `${sampleName} (NTC)` : sampleName,
            color: sampleColor(si),
          };
          bySample.set(si, trace);
        }
        trace.x.push(data.x[i]);
        trace.y.push(y);
        trace.text.push(hoverText(data, i, group, metric));
      }

      for (const trace of bySample.values()) {
        const order = trace.x.map((_, i) => i);
        if (sortMode === 'y_desc') order.sort((a, b) => (trace.y[b] ?? -Infinity) - (trace.y[a] ?? -Infinity));
        else if (sortMode === 'y_asc') order.sort((a, b) => (trace.y[a] ?? Infinity) - (trace.y[b] ?? Infinity));
        else order.sort((a, b) => trace.x[a] - trace.x[b]);
        trace.x = order.map(i => trace.x[i]);
        trace.y = order.map(i => trace.y[i]);
        trace.text = order.map(i => trace.text[i]);
      }

      let plotted = 0;
      for (const t of bySample.values()) plotted += t.x.length;
      // NTC traces draw a connecting line, which scattergl renders differently;
      // keep the whole figure on SVG unless it is genuinely large.
      const traceType = (HAS_WEBGL && plotted > GL_POINT_THRESHOLD) ? 'scattergl' : 'scatter';

      const traces = [...bySample.values()].sort((a, b) => a.si - b.si).map(t => {
        const styleAsNtc = t.isNtc && toggleNtc.checked;
        return {
          type: traceType,
          mode: styleAsNtc ? 'lines+markers' : (meta.samples.length <= 20 ? 'lines+markers' : 'markers'),
          name: t.name,
          x: t.x, y: t.y, text: t.text,
          hovertemplate: '%{text}<extra></extra>',
          marker: {
            size: styleAsNtc ? 6 : 7,
            opacity: styleAsNtc ? 0.55 : 0.86,
            color: t.color,
            symbol: styleAsNtc ? 'x' : sampleSymbolByIndex.get(t.si),
          },
          line: styleAsNtc ? { color: t.color, width: 1.2, dash: 'dash' } : (meta.samples.length <= 20 ? { color: t.color, width: 1.5 } : undefined),
          // No legendgroup: Plotly keeps same-group entries in one column, which
          // turns a 30-sample horizontal legend into a tall scrolling strip.
          // groupclick was already 'toggleitem', so the group bought nothing.
          legendrank: styleAsNtc ? 9000 : 1000,
        };
      });

      const reverseAxis = toggleStrand.checked && data.strand === '-';
      const xaxis = {
        title: meta.has_gene_data ? 'Exon Number' : 'Target Order',
        showgrid: false, zeroline: false, automargin: true,
        autorange: reverseAxis ? 'reversed' : true,
      };
      if (meta.has_gene_data && data.tv && data.tv.length) {
        xaxis.tickmode = 'array';
        xaxis.tickvals = data.tv;
        xaxis.ticktext = data.tt;
        xaxis.tickangle = data.tt.length > 18 ? -90 : 0;
        xaxis.tickfont = { size: data.tt.length > 30 ? 10 : 11 };
      }

      const strandLabel = data.strand === '-' ? 'strand: - (reversed)' : data.strand === '+' ? 'strand: +' : '';
      if (strandLabel) {
        strandChip.hidden = false;
        strandChip.textContent = strandLabel;
        strandChip.className = `chip ${data.strand === '-' ? 'strand-minus' : 'strand-plus'}`;
      } else {
        strandChip.hidden = true;
      }

      // A vertical legend clips once there are many samples, so spill it below the plot.
      const tc = themePlotColors();
      const manySamples = meta.samples.length > 12;
      const layout = {
        title: { text: `${group}${data.strand ? `  (${data.strand} strand)` : ''}`, font: { color: tc.fontColor, size: 17 } },
        font: { family: '"Segoe UI Variable","Aptos","SF Pro Display","Segoe UI",sans-serif', color: tc.fontColor },
        paper_bgcolor: 'rgba(0,0,0,0)',
        plot_bgcolor: tc.plotBg,
        hovermode: 'closest',
        hoverlabel: { bgcolor: tc.hoverBg, bordercolor: tc.hoverBorder, font: { color: tc.hoverFont, size: 12 } },
        showlegend: true,
        legend: manySamples
          ? { orientation: 'h', y: -0.28, yanchor: 'top', x: 0, font: { size: 10, color: tc.fontColor }, groupclick: 'toggleitem' }
          : { orientation: 'v', x: 1.09, xanchor: 'left', font: { color: tc.fontColor }, groupclick: 'toggleitem' },
        margin: { l: 70, r: manySamples ? 95 : 110, t: 60,
                  b: (meta.has_gene_data ? 120 : 70) + (manySamples ? 40 + Math.ceil(meta.samples.length / 6) * 16 : 0) },
        xaxis,
        yaxis: { title: metricLabel(metric), gridcolor: tc.gridColor,
                 zeroline: metric !== 'normalized_ratio', zerolinecolor: tc.zeroLineColor },
        shapes: buildShapes(data, metric),
        annotations: buildAnnotations(metric),
      };

      // Plotly.react does not reliably swap an SVG plot for a WebGL one in place.
      if (window.__lastTraceType && window.__lastTraceType !== traceType) Plotly.purge('plot');
      window.__lastTraceType = traceType;
      Plotly.react('plot', traces, layout, { responsive: true, displaylogo: false });
      renderCalls(group, data);
    }

    function currentGroup() {
      return (groupInput.value || '').trim() || meta.default_group || meta.groups[0] || '';
    }

    function resetControls() {
      groupInput.value = meta.default_group || meta.groups[0] || '';
      metricSelect.value = 'normalized_ratio';
      sortSelect.value = 'x';
      toggleRefs.checked = true;
      toggleBands.checked = true;
      toggleStrand.checked = true;
      toggleNtc.checked = true;
      render(currentGroup());
    }

    renderBtn.addEventListener('click', () => render(currentGroup()));
    resetBtn.addEventListener('click', resetControls);
    groupInput.addEventListener('keydown', e => { if (e.key === 'Enter') render(currentGroup()); });
    groupInput.addEventListener('input', () => { const v = (groupInput.value || '').trim(); if (groupIndex.has(v)) render(v); });
    for (const el of [toggleRefs, toggleBands, toggleStrand, toggleNtc, metricSelect, sortSelect]) {
      el.addEventListener('change', () => render(currentGroup()));
    }

    resetControls();
    window.__cnvReady = true;
