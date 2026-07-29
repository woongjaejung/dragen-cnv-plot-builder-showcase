/* Browser runtime for the synthetic per-gene CNV plot. */
(function () {
  "use strict";

  var config = JSON.parse(document.getElementById("cnv-config").textContent);
  var scriptByGene = new Map();
  var payloadCache = new Map();
  Object.keys(config.groupIds).forEach(function (gene) {
    scriptByGene.set(gene, document.getElementById(config.groupIds[gene]));
  });
  var geneNames = Object.keys(config.groupIds);

  var dom = {
    search: document.getElementById("gene-search"),
    metric: document.getElementById("metric"),
    plot: document.getElementById("plot"),
    status: document.getElementById("status"),
    currentGene: document.getElementById("current-gene"),
    strandChip: document.getElementById("strand-chip"),
    callsTable: document.getElementById("calls-table"),
    definitions: document.getElementById("definitions"),
    reference: document.getElementById("toggle-reference"),
    bands: document.getElementById("toggle-bands"),
    strand: document.getElementById("toggle-strand"),
    ntc: document.getElementById("toggle-ntc")
  };

  var state = {
    gene: config.defaultGroup,
    metric: "normalized_ratio",
    reference: true,
    bands: true,
    strand: true,
    ntc: true
  };

  var metricLabels = {
    normalized_ratio: "normalized ratio",
    tn_log2: "tn_log2",
    segment_mean: "segment mean"
  };

  var referenceValues = {
    normalized_ratio: [
      { value: 1.0, label: "diploid 1.0", color: "#53636f" },
      { value: 0.5, label: "het del 0.5", color: "#b55b62" },
      { value: 1.5, label: "het dup 1.5", color: "#397f8c" }
    ],
    tn_log2: [
      { value: 0.0, label: "diploid 0.0", color: "#53636f" },
      { value: -1.0, label: "het del -1.0", color: "#b55b62" },
      { value: 0.585, label: "het dup +0.585", color: "#397f8c" }
    ],
    segment_mean: []
  };

  var normalPalette = [
    "#176b87", "#d36c3f", "#567a4c", "#8d5a9e", "#b88328",
    "#247d78", "#a94b62", "#596c9e", "#8a6a42", "#3c6472"
  ];

  var filterDefinitions = {
    PASS: "caller accepted the call",
    cnvQual: "call quality below threshold",
    cnvLength: "segment shorter than the minimum reportable length",
    sampleQual: "sample-level QC below threshold"
  };

  function getPayload(gene) {
    if (!payloadCache.has(gene)) {
      var script = scriptByGene.get(gene);
      if (!script) {
        return null;
      }
      /* Parsing happens only when a gene is first viewed. */
      payloadCache.set(gene, JSON.parse(script.textContent));
    }
    return payloadCache.get(gene);
  }

  function escapeHtml(value) {
    return String(value)
      .replace(/&/g, "&amp;")
      .replace(/</g, "&lt;")
      .replace(/>/g, "&gt;")
      .replace(/"/g, "&quot;")
      .replace(/'/g, "&#039;");
  }

  function formatNumber(value) {
    if (value === null || value === undefined || Number.isNaN(value)) {
      return "n/a";
    }
    return Number(value).toFixed(3).replace(/\.000$/, "");
  }

  function formatMetricValue(value) {
    return formatNumber(value);
  }

  function displaySample(sample) {
    return sample.ntc ? sample.sample_id + " [NTC]" : sample.sample_id;
  }

  function orderedSamples(payload) {
    return payload.samples.slice().sort(function (left, right) {
      if (left.ntc !== right.ntc) {
        return left.ntc ? 1 : -1;
      }
      return left.sample_id.localeCompare(right.sample_id);
    });
  }

  function sampleColorMap(samples) {
    var colors = new Map();
    var normalIndex = 0;
    var ntcIndex = 0;
    samples.forEach(function (sample) {
      if (sample.ntc) {
        colors.set(sample.sample_id, ntcIndex === 0 ? "#81909b" : "#aab5bd");
        ntcIndex += 1;
      } else {
        colors.set(sample.sample_id, normalColor(normalIndex));
        normalIndex += 1;
      }
    });
    return colors;
  }

  function normalColor(index) {
    if (index < normalPalette.length) {
      return normalPalette[index];
    }
    var hue = (index * 137.508) % 360;
    return "hsl(" + hue.toFixed(1) + ", 55%, 42%)";
  }

  function targetOrder(payload) {
    var indices = payload.targets.x.map(function (_, index) { return index; });
    if (state.strand && payload.strand === "-") {
      indices.reverse();
    }
    return indices;
  }

  function metricValue(sample, index) {
    if (state.metric === "normalized_ratio") {
      return Math.pow(2, sample.tn_log2[index]);
    }
    if (state.metric === "tn_log2") {
      return sample.tn_log2[index];
    }
    return sample.segment_mean[index];
  }

  function callsForTarget(payload, sampleId, index) {
    var start = payload.targets.start[index];
    var end = payload.targets.end[index];
    return payload.calls.filter(function (call) {
      return call.sample_id === sampleId &&
        start <= call.end && end >= call.start;
    });
  }

  function callHoverText(calls) {
    if (!calls.length) {
      return "none";
    }
    return calls.map(function (call) {
      return call.type + " CN=" + call.copy_number + " " + call.filter +
        " QUAL=" + formatNumber(call.qual);
    }).join("; ");
  }

  function targetHoverData(payload, sample, order) {
    return order.map(function (index) {
      var value = metricValue(sample, index);
      var calls = callsForTarget(payload, sample.sample_id, index);
      return [
        displaySample(sample),
        payload.gene,
        payload.targets.exon_number[index],
        payload.targets.target_name[index],
        payload.contig + ":" + payload.targets.start[index] + "-" + payload.targets.end[index],
        metricLabels[state.metric] + " = " + formatMetricValue(value),
        callHoverText(calls)
      ];
    });
  }

  function boundaryShapes(payload, order) {
    var shapes = [];
    for (var position = 1; position < order.length; position += 1) {
      var previous = payload.targets.exon_number[order[position - 1]];
      var current = payload.targets.exon_number[order[position]];
      if (previous !== current) {
        shapes.push({
          type: "line",
          xref: "x",
          yref: "paper",
          x0: position - 0.5,
          x1: position - 0.5,
          y0: 0,
          y1: 1,
          line: { color: "#d6dee3", width: 1, dash: "dot" },
          layer: "below"
        });
      }
    }
    return shapes;
  }

  function referenceShapes() {
    return referenceValues[state.metric].map(function (reference) {
      return {
        type: "line",
        xref: "paper",
        yref: "y",
        x0: 0,
        x1: 1,
        y0: reference.value,
        y1: reference.value,
        line: { color: reference.color, width: 1, dash: "dash" },
        layer: "below"
      };
    });
  }

  function referenceAnnotations() {
    return referenceValues[state.metric].map(function (reference) {
      return {
        xref: "paper",
        yref: "y",
        x: 0,
        y: reference.value,
        xanchor: "left",
        yanchor: "bottom",
        text: reference.label,
        showarrow: false,
        font: { size: 10, color: reference.color },
        bgcolor: "rgba(255,255,255,0.72)",
        borderpad: 2
      };
    });
  }

  function passBandShapes(payload, order) {
    var seen = new Set();
    var shapes = [];
    payload.calls.forEach(function (call) {
      if (call.filter !== "PASS") {
        return;
      }
      var visible = order.filter(function (index) {
        return payload.targets.start[index] <= call.end &&
          payload.targets.end[index] >= call.start;
      });
      if (!visible.length) {
        return;
      }
      var key = String(call.start) + ":" + String(call.end);
      if (seen.has(key)) {
        return;
      }
      seen.add(key);
      var positions = visible.map(function (index) { return order.indexOf(index); });
      var isDup = call.type === "DUP";
      shapes.push({
        type: "rect",
        xref: "x",
        yref: "paper",
        x0: Math.min.apply(null, positions) - 0.5,
        x1: Math.max.apply(null, positions) + 0.5,
        y0: 0,
        y1: 1,
        fillcolor: isDup ? "rgba(57,127,140,0.08)" : "rgba(181,91,98,0.08)",
        line: { width: 0 },
        layer: "below"
      });
    });
    return shapes;
  }

  function exonTicks(payload, order) {
    var positionsByExon = new Map();
    order.forEach(function (index, position) {
      var exon = payload.targets.exon_number[index];
      if (!positionsByExon.has(exon)) {
        positionsByExon.set(exon, []);
      }
      positionsByExon.get(exon).push(position);
    });
    var values = [];
    var text = [];
    positionsByExon.forEach(function (positions, exon) {
      values.push(positions.reduce(function (total, position) {
        return total + position;
      }, 0) / positions.length);
      text.push(String(exon));
    });
    return { values: values, text: text };
  }

  function makeTraces(payload, order, samples, colors) {
    return samples.map(function (sample, sampleIndex) {
      var color = colors.get(sample.sample_id);
      var styledNtc = state.ntc && sample.ntc;
      var values = order.map(function (index) { return metricValue(sample, index); });
      return {
        x: order.map(function (_, index) { return index; }),
        y: values,
        type: "scatter",
        mode: "lines+markers",
        name: displaySample(sample),
        legendrank: sampleIndex,
        connectgaps: false,
        line: {
          color: color,
          width: styledNtc ? 1.5 : 2,
          dash: styledNtc ? "dash" : "solid"
        },
        marker: {
          color: color,
          size: styledNtc ? 8 : 6,
          symbol: styledNtc ? "x" : "circle",
          line: { width: styledNtc ? 1.2 : 0, color: color }
        },
        opacity: styledNtc ? 0.48 : 0.9,
        customdata: targetHoverData(payload, sample, order),
        hovertemplate:
          "<b>%{customdata[0]}</b><br>gene: %{customdata[1]}<br>" +
          "exon: %{customdata[2]}<br>target: %{customdata[3]}<br>" +
          "region: %{customdata[4]}<br>%{customdata[5]}<br>" +
          "overlapping calls: %{customdata[6]}<extra></extra>"
      };
    });
  }

  function plotLayout(payload, order) {
    var ticks = exonTicks(payload, order);
    var shapes = boundaryShapes(payload, order);
    if (state.reference) {
      shapes = shapes.concat(referenceShapes());
    }
    if (state.bands) {
      shapes = shapes.concat(passBandShapes(payload, order));
    }
    return {
      height: 520,
      margin: { l: 64, r: 220, t: 20, b: 70 },
      paper_bgcolor: "#ffffff",
      plot_bgcolor: "#ffffff",
      hovermode: "closest",
      hoverlabel: { bgcolor: "#17212b", font: { color: "#ffffff", size: 12 } },
      legend: {
        orientation: "v",
        x: 1.01,
        xanchor: "left",
        y: 1,
        yanchor: "top",
        font: { size: 11 },
        bgcolor: "rgba(255,255,255,0.88)",
        bordercolor: "#d9e1e7",
        borderwidth: 1
      },
      xaxis: {
        title: { text: "exon number (transcript order)" },
        tickmode: "array",
        tickvals: ticks.values,
        ticktext: ticks.text,
        range: [-0.5, order.length - 0.5],
        showgrid: false,
        zeroline: false,
        fixedrange: false
      },
      yaxis: {
        title: { text: metricLabels[state.metric] },
        zeroline: false,
        gridcolor: "#e9eef1",
        automargin: true
      },
      shapes: shapes,
      annotations: state.reference ? referenceAnnotations() : [],
      uirevision: payload.gene + "-" + state.metric
    };
  }

  function renderPlot(payload) {
    var order = targetOrder(payload);
    var samples = orderedSamples(payload);
    var colors = sampleColorMap(samples);
    var traces = makeTraces(payload, order, samples, colors);
    var layout = plotLayout(payload, order);
    Plotly.react(dom.plot, traces, layout, {
      responsive: true,
      displaylogo: false,
      modeBarButtonsToRemove: ["lasso2d", "select2d"]
    });
    dom.status.textContent = payload.samples.length + " samples overlaid; " +
      payload.targets.start.length + " capture targets across " +
      new Set(payload.targets.exon_number).size + " exons. " +
      (payload.strand === "-" && state.strand ? "Negative strand reversed 5' -> 3'." : "Genomic coordinate order shown.");
    return colors;
  }

  function exonRange(exons) {
    if (!exons || !exons.length) {
      return "n/a";
    }
    var sorted = exons.slice().sort(function (left, right) { return left - right; });
    var ranges = [];
    var start = sorted[0];
    var previous = sorted[0];
    for (var index = 1; index < sorted.length; index += 1) {
      if (sorted[index] !== previous + 1) {
        ranges.push(start === previous ? String(start) : start + "-" + previous);
        start = sorted[index];
      }
      previous = sorted[index];
    }
    ranges.push(start === previous ? String(start) : start + "-" + previous);
    return ranges.join(", ");
  }

  function renderCalls(payload, colors) {
    var calls = payload.calls.slice();
    if (!calls.length) {
      dom.callsTable.innerHTML = '<p class="empty">No calls overlap targets in this gene.</p>';
      dom.definitions.innerHTML = "";
      return;
    }
    var rows = calls.map(function (call) {
      var color = colors.get(call.sample_id) || "#81909b";
      var filterClass = call.filter === "PASS" ? "filter-tag pass" : "filter-tag";
      return "<tr>" +
        "<td><span class=\"sample-cell\"><span class=\"sample-swatch\" style=\"background:" + color + "\"></span>" +
        escapeHtml(call.sample_id) + "</span></td>" +
        "<td>" + escapeHtml(call.type) + "</td>" +
        "<td>" + escapeHtml(call.copy_number) + "</td>" +
        "<td><span class=\"" + filterClass + "\">" + escapeHtml(call.filter) + "</span></td>" +
        "<td>" + formatNumber(call.qual) + "</td>" +
        "<td>" + escapeHtml(payload.contig + ":" + call.start + "-" + call.end) + "</td>" +
        "<td>" + escapeHtml(exonRange(call.exons)) + "</td>" +
        "</tr>";
    }).join("");
    dom.callsTable.innerHTML =
      "<table><thead><tr><th>sample</th><th>type</th><th>copy number</th>" +
      "<th>filter tag</th><th>QUAL</th><th>call region</th>" +
      "<th>exons visible</th></tr></thead><tbody>" + rows + "</tbody></table>";

    var presentFilters = [];
    calls.forEach(function (call) {
      if (presentFilters.indexOf(call.filter) === -1) {
        presentFilters.push(call.filter);
      }
    });
    dom.definitions.innerHTML = presentFilters.map(function (filter) {
      var definition = filterDefinitions[filter] || "caller filter tag present in this view";
      return "<p><strong>" + escapeHtml(filter) + "</strong> = " + escapeHtml(definition) + "</p>";
    }).join("");
  }

  function updateHeader(payload) {
    dom.currentGene.textContent = payload.gene;
    var reversed = payload.strand === "-" && state.strand;
    dom.strandChip.innerHTML =
      '<span class="arrow">' + (payload.strand === "-" ? "-" : "+") + "</span> " +
      escapeHtml(payload.strand === "-" ? (reversed ? "negative / 5' -> 3'" : "negative / genomic") : "positive / genomic");
  }

  function render() {
    var payload = getPayload(state.gene);
    if (!payload) {
      dom.status.innerHTML = '<span class="error">Gene not found. Choose one from the search list.</span>';
      return;
    }
    updateHeader(payload);
    var colors = renderPlot(payload);
    renderCalls(payload, colors);
  }

  function chooseGene(value) {
    var normalized = value.trim().toLowerCase();
    var exact = geneNames.find(function (gene) {
      return gene.toLowerCase() === normalized;
    });
    if (exact) {
      state.gene = exact;
      dom.search.value = state.gene;
      render();
      return;
    }
    var prefix = geneNames.find(function (gene) {
      return gene.toLowerCase().indexOf(normalized) === 0;
    });
    if (prefix && normalized.length > 1) {
      state.gene = prefix;
      render();
    }
  }

  function bindControls() {
    dom.search.addEventListener("input", function () { chooseGene(dom.search.value); });
    dom.metric.addEventListener("change", function () {
      state.metric = dom.metric.value;
      render();
    });
    dom.reference.addEventListener("change", function () {
      state.reference = dom.reference.checked;
      render();
    });
    dom.bands.addEventListener("change", function () {
      state.bands = dom.bands.checked;
      render();
    });
    dom.strand.addEventListener("change", function () {
      state.strand = dom.strand.checked;
      render();
    });
    dom.ntc.addEventListener("change", function () {
      state.ntc = dom.ntc.checked;
      render();
    });
  }

  function initialize() {
    bindControls();
    dom.search.value = state.gene;
    if (typeof Plotly === "undefined") {
      dom.status.innerHTML = '<span class="error">Plotly did not load from the local vendor path.</span>';
      return;
    }
    render();
  }

  if (document.readyState === "loading") {
    document.addEventListener("DOMContentLoaded", initialize);
  } else {
    initialize();
  }
})();
