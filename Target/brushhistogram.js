// Brushing based histogram display


// Adapted from Mike Bostock's example: http://bl.ocks.org/mbostock/3048450, http://tributary.io/inlet/5998335 and http://embed.plnkr.co/agRZx6/script.js
var plotHist = (function(){

  var margin = {top: 25, right: 15, bottom: 30, left: 50},
      width = 400 - margin.left - margin.right,
      height = 300 - margin.top - margin.bottom;

  var svg = d3.select(".stat").append("div")
      .classed("histogram", true);

  var svg = svg
      .append("svg")
      .attr("width", width + margin.left + margin.right)
      .attr("height", height + margin.top + margin.bottom)
    .append("g")
      .attr("transform", "translate(" + margin.left + "," + margin.top + ")");

  var svgtitle = svg.append("text")
     .attr("x", (width / 2))
     .attr("y", 0 - (margin.top / 2))
     .attr("text-anchor", "middle");

  svg.append("g")
      .classed("x axis", true);

  svg.append("g")
      .classed("y axis", true);

  return function(values){
    if (values.length < 3) {
       return
    }

    svgtitle
     .text("Chiral Center");

    var x = d3.scale.linear()
        .domain([0, d3.max(values)])
        .range([0, width]);

    // Generate a histogram using twenty uniformly-spaced bins.
    var data = d3.layout.histogram()
        .bins(x.ticks(20))
        (values);

    var y = d3.scale.linear()
        .domain([0, d3.max(data, function(d) { return d.y; })])
        .range([height, 0]);

    var xAxis = d3.svg.axis()
        .scale(x)
        .orient("bottom");

    var yAxis = d3.svg.axis()
        .scale(y)
        .orient("left");

    var bar = svg.selectAll(".bar")
        .data(data);

    bar.exit().remove();

    bar.enter()
        .append("g")
          .attr("class", "bar")
          .attr("transform", function(d) { return "translate(" + x(d.x) + "," + y(d.y) + ")"; })
        .append("rect")
          .attr("x", 1)
          .attr("width", x(data[0].dx) - 1)
          .attr("height", function(d) { return height - y(d.y); });

    bar.attr("transform", function(d) { return "translate(" + x(d.x) + "," + y(d.y) + ")"; })
        .select("rect")
         .attr("width", x(data[0].dx) - 1)
         .attr("height", function(d) { return height - y(d.y); });

    svg.selectAll("g.x.axis")
        .attr("transform", "translate(0," + height + ")")
        .call(xAxis);

    svg.selectAll('g.y.axis')
        .call(yAxis);
  }
})()


var brushing = function(){
  var brush = graph.append("g")
                .attr("class", "brush");

  brush
    .call(d3.svg.brush()
      .x(xrange)
      .y(yrange)
      .on("brush", function() {
        var extent = d3.event.target.extent();
        var values = [];
        node.classed("selectedbrush", function(d) {
          if ( extent[0][0] <= d.x && d.x < extent[1][0]
               && extent[0][1] <= d.y && d.y < extent[1][1] && !d.children) {
              values.push(d.continuous.chiral);
              return true;
           }
          });
       // plotHist(values);
      }));
}

// Adapted from Crossfilter example
var chart;
var charts = [];
var nodes;
var dimen_list = [];
var group_list = [];
var plotHists = function(){
    nodes  = flatten(root);
    nodes = nodes.filter(function(d){ return d.children == null } );
    var contin = crossfilter(nodes); // Continous data
    // Dimensions and groups
    var keys = Object.keys(nodes[0].continuous);
    console.log(keys);
    for (var i = 0; i < keys.length; i++) {
        dimen_list.push(contin.dimension(function(d){ return d.continuous[keys[i]]}));
    }
    for (var i = 0; i < dimen_list.length; i++) {
        group_list.push(dimen_list[i].group(Math.floor));
    }

    // Charts
    for (var i in d3.range(Object.keys(nodes[0].continuous).length)) {
        console.log(i);
        console.log([dimen_list[i].bottom(1)[0].continuous[keys[i]], dimen_list[i].top(1)[0].continuous[keys[i]]]);
        charts.push(barChart()
                        .dimension(dimen_list[i])
                        .group(group_list[i])
                      .x(d3.scale.linear()
                        .domain([dimen_list[i].bottom(1)[0].continuous[keys[i]], dimen_list[i].top(1)[0].continuous[keys[i]]])
                        .rangeRound([0, 10*24])))
    }
    // listen to the chart's brush events to update the display.
    var render = function(method){
        d3.select(this).call(method);
    }

    // re-render everything
    var renderAll = function(){
        chart.each(render);
    }


    var svg = d3.select(".stat").append("div")
        .classed("histogram", true);

    chart = svg.selectAll(".chart")
        .data(charts)
        .each(function(chart) { chart.on("brush", renderAll).on("brushend", renderAll); });

    chart.enter();

    renderAll();

    function barChart() {
      if (!barChart.id) barChart.id = 0;

      var margin = {top: 10, right: 10, bottom: 20, left: 10},
          x,
          y = d3.scale.linear().range([100, 0]),
          id = barChart.id++,
          axis = d3.svg.axis().orient("bottom"),
          brush = d3.svg.brush(),
          brushDirty,
          dimension,
          group,
          round;

      function chart(div) {
        var width = x.range()[1],
            height = y.range()[0];

        y.domain([0, group.top(1)[0].value]);

        div.each(function() {
          var div = d3.select(this),
              g = div.select("g");

          // Create the skeletal chart.
          if (g.empty()) {
            div.select(".title").append("a")
                .attr("href", "javascript:reset(" + id + ")")
                .attr("class", "reset")
                .text("reset")
                .style("display", "none");

            g = div.append("svg")
                .attr("width", width + margin.left + margin.right)
                .attr("height", height + margin.top + margin.bottom)
              .append("g")
                .attr("transform", "translate(" + margin.left + "," + margin.top + ")");

            g.append("clipPath")
                .attr("id", "clip-" + id)
              .append("rect")
                .attr("width", width)
                .attr("height", height);

            g.selectAll(".bar")
                .data(["background", "foreground"])
              .enter().append("path")
                .attr("class", function(d) { return d + " bar"; })
                .datum(group.all());

            g.selectAll(".foreground.bar")
                .attr("clip-path", "url(#clip-" + id + ")");

            g.append("g")
                .attr("class", "axis")
                .attr("transform", "translate(0," + height + ")")
                .call(axis);

            // Initialize the brush component with pretty resize handles.
            var gBrush = g.append("g").attr("class", "brush").call(brush);
            gBrush.selectAll("rect").attr("height", height);
            gBrush.selectAll(".resize").append("path").attr("d", resizePath);
          }

          // Only redraw the brush if set externally.
          if (brushDirty) {
            brushDirty = false;
            g.selectAll(".brush").call(brush);
            div.select(".title a").style("display", brush.empty() ? "none" : null);
            if (brush.empty()) {
              g.selectAll("#clip-" + id + " rect")
                  .attr("x", 0)
                  .attr("width", width);
            } else {
              var extent = brush.extent();
              g.selectAll("#clip-" + id + " rect")
                  .attr("x", x(extent[0]))
                  .attr("width", x(extent[1]) - x(extent[0]));
            }
          }

          g.selectAll(".bar").attr("d", barPath);
        });

        function barPath(groups) {
          var path = [],
              i = -1,
              n = groups.length,
              d;
          while (++i < n) {
            d = groups[i];
            path.push("M", x(d.key), ",", height, "V", y(d.value), "h9V", height);
          }
          return path.join("");
        }

        function resizePath(d) {
          var e = +(d == "e"),
              x = e ? 1 : -1,
              y = height / 3;
          return "M" + (.5 * x) + "," + y
              + "A6,6 0 0 " + e + " " + (6.5 * x) + "," + (y + 6)
              + "V" + (2 * y - 6)
              + "A6,6 0 0 " + e + " " + (.5 * x) + "," + (2 * y)
              + "Z"
              + "M" + (2.5 * x) + "," + (y + 8)
              + "V" + (2 * y - 8)
              + "M" + (4.5 * x) + "," + (y + 8)
              + "V" + (2 * y - 8);
        }
      }

      brush.on("brushstart.chart", function() {
        var div = d3.select(this.parentNode.parentNode.parentNode);
        div.select(".title a").style("display", null);
      });

      brush.on("brush.chart", function() {
        var g = d3.select(this.parentNode),
            extent = brush.extent();
        if (round) g.select(".brush")
            .call(brush.extent(extent = extent.map(round)))
          .selectAll(".resize")
            .style("display", null);
        g.select("#clip-" + id + " rect")
            .attr("x", x(extent[0]))
            .attr("width", x(extent[1]) - x(extent[0]));
        dimension.filterRange(extent);
      });

      brush.on("brushend.chart", function() {
        if (brush.empty()) {
          var div = d3.select(this.parentNode.parentNode.parentNode);
          div.select(".title a").style("display", "none");
          div.select("#clip-" + id + " rect").attr("x", null).attr("width", "100%");
          dimension.filterAll();
        }
      });

      chart.margin = function(_) {
        if (!arguments.length) return margin;
        margin = _;
        return chart;
      };

      chart.x = function(_) {
        if (!arguments.length) return x;
        x = _;
        axis.scale(x);
        brush.x(x);
        return chart;
      };

      chart.y = function(_) {
        if (!arguments.length) return y;
        y = _;
        return chart;
      };

      chart.dimension = function(_) {
        if (!arguments.length) return dimension;
        dimension = _;
        return chart;
      };

      chart.filter = function(_) {
        if (_) {
          brush.extent(_);
          dimension.filterRange(_);
        } else {
          brush.clear();
          dimension.filterAll();
        }
        brushDirty = true;
        return chart;
      };

      chart.group = function(_) {
        if (!arguments.length) return group;
        group = _;
        return chart;
      };

      chart.round = function(_) {
        if (!arguments.length) return round;
        round = _;
        return chart;
      };

      return d3.rebind(chart, brush, "on");
    };
};
