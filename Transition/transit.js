// custom d3 javascript

// data comes from uploaded data to vida.io
// if you have your own data, overwrite its value
// drawing area is #canvas, define var svg = d3.select("#canvas")

var width = 900,
    height = 800;

var color = d3.scale.category20();

var force0 = d3.layout.force()
    .charge(-120)
    .linkDistance(30)
    .size([width, height]);

var force1 = d3.layout.force()
    .charge(-120)
    .linkDistance(30)
    .size([width, height]);

var svg = d3.select("#canvas").append("svg")
    .attr("width", width)
    .attr("height", height);

var link, node, nodes0, nodes1, links0, links1, linkG = svg.append("g");

var drawD3Document = function(data) {
  nodes0 = flatten(data[0]),
  links0 = d3.layout.tree().links(nodes0),
  nodes1 = flatten(data[1]),
  links1 = d3.layout.tree().links(nodes1);

  force0
      .nodes(nodes0)
      .links(links0)
      .start();

  force1
      .nodes(nodes1)
      .links(links1)
      .start();

  link = linkG.selectAll(".link")
      .data(links0, function(d) { return d.target.id; })
    .enter().append("line")
      .attr("class", "link");

  node = svg.selectAll(".node")
      .data(nodes0)
    .enter().append("circle")
      .attr("class", "node")
      .attr("r", 5)
      .style("fill", function(d) { return color(d.group); });

  node.append("title")
      .text(function(d) { return d.name; });

  force0.on("tick", function() {
    link.attr("x1", function(d) { return d.source.x; })
        .attr("y1", function(d) { return d.source.y; })
        .attr("x2", function(d) { return d.target.x; })
        .attr("y2", function(d) { return d.target.y; });

    node.attr("cx", function(d) { return d.x; })
        .attr("cy", function(d) { return d.y; });
  });

  setTimeout(function() {
    force0.stop();
    link.transition().duration(1000).style("stroke-opacity", 0).remove();

    node = svg.selectAll(".node")
      .data(nodes1, function(d) { return d.name; });
    node.enter().append("circle")
      .attr("class", "node")
      .attr("r", 5)
      .style("fill", function(d) { return color(d.group); });
    node.exit().remove();

    node.transition().delay(1000).duration(1000)
        .attr("cx", function(d) { return d.x; })
        .attr("cy", function(d) { return d.y; });

    setTimeout(function() {
      link = linkG.selectAll(".link")
        .data(links1)
      .enter().append("line")
        .attr("class", "link")
        .style("stroke-width", function(d) { return Math.sqrt(d.value); })
        .style("stroke-opacity", 0);
      link
        .transition().duration(1000)
        .style("stroke-opacity", 0.6)
        .attr("x1", function(d) { return d.source.x; })
        .attr("y1", function(d) { return d.source.y; })
        .attr("x2", function(d) { return d.target.x; })
        .attr("y2", function(d) { return d.target.y; });

      force1.on("tick", function() {
        link.attr("x1", function(d) { return d.source.x; })
            .attr("y1", function(d) { return d.source.y; })
           .attr("x2", function(d) { return d.target.x; })
           .attr("y2", function(d) { return d.target.y; });

       node.attr("cx", function(d) { return d.x; })
           .attr("cy", function(d) { return d.y; });
     });
    }, 2000);
  }, 5000);

};

function flatten(root) {
  var nodes = [], i = 0;

  function recurse(node) {
    if (node.children) node.size = node.children.reduce(function(p, v) { return p + recurse(v); }, 0);
    if (!node.id) node.id = ++i;
    nodes.push(node);
    return node.size;
  }

  root.size = recurse(root);
  return nodes;
}
