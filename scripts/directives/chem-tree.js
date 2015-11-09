/**
 * Created by ajing on 9/11/15.
 */

'use strict';

/**
 * @ngdoc directive
 * @name frontendApp.directive:chem-tree
 * @description
 * # chem-tree
 */

angular.module('frontendApp')
  .directive('chemTree', function () {
      function link(scope, elements) {
          var xScale, yScale, activityScale, activityColor, slogpScale, ligeffScale,
            borderScale, borderColor, sizeScale, force, linkDOM, nodeDOM;

          //setup
          var el = elements[0];

          //append the svg element
          var svg = d3.select(el)
            .append('svg')
            .attr({class: 'viz'});

          //clicking anywhere should set selected to none.  This should default away if clicking on an object
          svg.on('click', function(){
              if (d3.event.defaultPrevented) { return; }
              scope.selected = null;
              scope.$apply();
          });

          var vis = svg.append('g'); // the zoom container

          d3.selectAll('.viz').append('g').append('svg')
            .attr('x', '60')
            .attr('y', '100')
            .append('g')
            .attr('transform', 'translate(0, 10)').classed('colorbarA',true);// color bar Activyty

          d3.selectAll('.viz').append('g').append('svg')
            .attr('x', '20')
            .attr('y', '100')
            .append('g')
            .attr('transform', 'translate(0, 10)').classed('colorbarB',true);// color bar border

          function flatten(root){
              var nodes = [], i = 0;

              function recurse(node) {
                  if (node.children) { node.size = node.children.reduce(function(p, v) { return p + recurse(v); }, 0); }
                  nodes.push(node);
                  return node.size;
              }

              root.size = recurse(root);
              return nodes;
          }

          function tick() {
              linkDOM
                .attr('x1', function(d) { return xScale(d.source.x); })
                .attr('y1', function(d) { return yScale(d.source.y); })
                .attr('x2', function(d) { return xScale(d.target.x); })
                .attr('y2', function(d) { return yScale(d.target.y); });

              nodeDOM
                .attr('cx', function(d) { return xScale(d.x); })
                .attr('cy', function(d) { return yScale(d.y); });
          }

         activityScale = d3.scale.linear()
            .domain([4, 9])
            .clamp(true)
            .range(['hsl(300,80%,50%)', 'hsl(0,80%,50%)'])
            .interpolate(d3.interpolateString);

        /*
            .range(['#008000', '#FFFF00', '#FF0000']);
            .interpolate(d3.interpolateRgb)
         */

          slogpScale = d3.scale.linear()
            .domain([5, -5])
            .clamp(true)
            .range(['hsl(300,80%,50%)', 'hsl(0,80%,50%)'])
            .interpolate(d3.interpolateString);

          ligeffScale = d3.scale.linear()
            .domain([0, 0.5])
            .clamp(true)
            .range(['hsl(300,80%,50%)', 'hsl(0,80%,50%)'])
            .interpolate(d3.interpolateString);


          //update data
          scope.$watch( 'treeType', function(newTreeType) {
            //function(scope) { return scope.treeType === null; }, function() {

              if (scope.treeType === undefined || scope.data === undefined) {
                  return;
              }

              console.log('test');
              var colorExtent, colorbarB; // for border color
              var click, update; // for clicking

              var root  = scope.data.trees[scope.treeType],
                nodes = flatten(root),
                links = d3.layout.tree().links(nodes);

              force = scope.data.forces[scope.treeType];

            //extract the scales
              xScale = d3.scale.linear()
                .domain(d3.extent(nodes, function(d) { return d.x; }))
                .range([0.1 * window.innerWidth, 0.9 * window.innerWidth]);

              yScale = d3.scale.linear()
                .domain(d3.extent(nodes, function(d) { return d.y; }))
                .range([0.1 * window.innerHeight, 0.9 * window.innerHeight]);

              sizeScale = d3.scale.linear()
                .domain(d3.extent(nodes, function(d) { if (d.name[0] === 'B') return d.size; }))
                .range([4, 10]);

              var colorbarA = colorBar()
                .color(activityScale).size(350).lineWidth(20).precision(4).tickFormat(d3.format('g'));

              d3.select('.colorbarA').call(colorbarA);


              if (scope.circleBorderType === 'SLogP') {
                borderScale = slogpScale;
                colorbarB = colorBar()
                  .color(borderScale).size(350).lineWidth(20).precision(4).tickFormat(d3.format('g'));
              } else if (scope.circleBorderType === 'Lig_Eff') {
                borderScale = ligeffScale;
                colorbarB = colorBar()
                  .color(borderScale).size(350).lineWidth(20).precision(4).tickFormat(d3.format('g'));
              } else {
                colorExtent = d3.extent(nodes, function(d) { if (d.name[0] === 'B') { return d.stroke; } });
                borderScale = d3.scale.linear()
                  .domain([colorExtent[0], d3.mean(colorExtent), colorExtent[1]])
                  .range(['#008000', '#FFFF00', '#FF0000']);
              }

              if (scope.circleBorderType === 'None') {
                d3.select('.colorbarB').style('visibility','hidden');
              } else {
                d3.select('.colorbarB').style('visibility','visible').call(colorbarB);
              }


              activityColor = function(d) {
                  return d._children ? '#3182bd' : d.children ? '#D3D3D3' : activityScale(d.fill);
              };

              borderColor = function(d) {
                //console.log(d.stroke);
                //if (d.stroke == 0) { return '#c6dbef'; }
                return d._children ? '#CCC' : d.children ? '#D3D3D3' : borderScale(d.stroke);
              };



              /** Drag behavior configuration **/
              // zoom and drag may interfere with each other, so here redefine drag function
              var drag = d3.behavior.drag()
                .origin(function(d){ return d; }) // identify function
                .on('dragstart', dragstarted)
                .on('drag', dragged)
                .on('dragend', dragended);

              function dragstarted(d){
                  d3.event.sourceEvent.stopPropagation();
                  d3.select(this).classed('dragged', true);
              }

              function dragged(d){
                  var mouselocation = d3.mouse(vis.node());

                  d.x = xScale.invert(mouselocation[0]);
                  d.y = yScale.invert(mouselocation[1]);
                  tick(); // re-position this node and links connected to this node
              }

              function dragended(d){
                  d3.select(this).classed('dragging', false);
                  force.resume();
              }


              /** Click behavior configuration **/
              // Toggle children on click.
              function clickOnInterNode(d) {
                if (d.children) {
                  d._children = d.children;
                  d.children = null;
                } else {
                  d.children = d._children;
                  d._children = null;
                }
                force.start();
                update();
              }

              function click(compound) {
                if (compound.name[0] !== 'B') { return clickOnInterNode(compound); }
                if (d3.event.defaultPrevented) { return; }
                d3.event.preventDefault();
                if (compound === scope.selected) {
                  scope.selected = null;
                } else {
                  scope.selected = scope.data.compounds[parseInt(compound.name.substring(1))];
                }
                force.start();
                scope.$apply();
              }



              function update() {
                var root  = scope.data.trees[scope.treeType],
                  nodes = flatten(root),
                  links = d3.layout.tree().links(nodes);

                //create the selections and bind them to the data
                //console.log(links);
                nodeDOM = vis.selectAll('circle').data(nodes, function(d) { return d.name; });


                // Enter any new nodes.
                nodeDOM.enter().append('svg:circle')
                  .attr('class', 'node')
                  .attr('cx', function(d) { return xScale(d.x); })
                  .attr('cy', function(d) { return yScale(d.y); })
                  .attr('r', function(d) { return d.children ? 2 : sizeScale(d.size); })
                  .style('fill', activityColor)
                  .style('stroke', borderColor)
                  .style('stroke-width', function(d) { return d.strokeWidth; })
                  .on('click', click)
                  // .on('mouseover', mouseover)
                  .call(drag) // attach drag behavior to new circles
                  .append('svg:title')
                  .text( function(d){ return d.name; });
                // .classed('selected', function (d) { return d.name === scope.id; });


                // transition from old to new
                nodeDOM
                  //.transition().duration(750)
                  .attr('cx', function(d) { return xScale(d.x); })
                  .attr('cy', function(d) { return yScale(d.y); })
                  .attr('r', function(d) { return d.children ? 2 : sizeScale(d.size); });
                //nodeDOM.filter(function(d) { return d.name[0] === 'B' ? null : this; })
                //  .attr('cx', function(d) { return xScale(d.x); })
                //  .attr('cy', function(d) { return yScale(d.y); })
                //  .attr('r', function(d) { return sizeScale(d.r); });



                // Exit any old nodes.
                nodeDOM.exit().remove();

                // 2. update new links
                linkDOM = vis.selectAll('line')
                  .data(links, function(d) { return d.target.name; });

                // transition from old to new
                linkDOM
                  //     .transition(0).duration(750)
                  .attr('x1', function(d) { return xScale(d.source.x); })
                  .attr('y1', function(d) { return yScale(d.source.y); })
                  .attr('x2', function(d) { return xScale(d.target.x); })
                  .attr('y2', function(d) { return yScale(d.target.y); });

                // Enter any new links.
                linkDOM.enter().insert('svg:line', '.node')
                  .attr('class', 'link')
                  .attr('x1', function(d) { return xScale(d.source.x); })
                  .attr('y1', function(d) { return yScale(d.source.y); })
                  .attr('x2', function(d) { return xScale(d.target.x); })
                  .attr('y2', function(d) { return yScale(d.target.y); })
                  .style('stroke', '#D3D3D3')
                  .style('stroke-width', '5px');

                // 2. Exit previous links
                linkDOM.exit().remove();

                // 3. transition
                //linkDOM.transition().duration(750)
                //  .style('stroke-opacity', 0.5);

              }

              update();


              //zoomer
              /** Zoom behavior configuration **/

              function zoom() {
                  tick(); // update position by tick, so the actual d.x, d.y won't change
              }

              var zoomer = d3.behavior.zoom()
                  // allow only 10 times zoom in or out
                .scaleExtent([0.1, 10])
                  // attach zoom function for variable modification
                .on('zoom', zoom);

              // let zoomer ajust coordinates by xScale and yScale
              zoomer.x(xScale).y(yScale);

              svg.call(zoomer);

              force.on('tick', tick);

              // start the force
              console.log(force.size());

              // set default gravity
              //scope.gravityValue = force.gravity;
          });

/*          scope.$watch('treeType', function(newTreeType, oldTreeType) {
              if (scope.treeType === null || scope.data === undefined || oldTreeType === null) {
                  return;
              }

              var root  = scope.data.trees[scope.treeType],
                nodes = flatten(root),
                links = d3.layout.tree().links(nodes);

              // stop the previous force
              if (oldTreeType !== null) {
                  scope.data.forces[oldTreeType].stop();
              }



              // transition from old to new
              nodeDOM.filter(function(d) { return d.name[0] === 'B' ? this : null; })
                .transition().delay(100).duration(750)
                .attr('cx', function(d) { return xScale(d.x); })
                .attr('cy', function(d) { return yScale(d.y); });
              nodeDOM.filter(function(d) { return d.name[0] === 'B' ? null : this; })
                .attr('cx', function(d) { return xScale(d.x); })
                .attr('cy', function(d) { return yScale(d.y); });
          });*/


          scope.$watch('circleSizeType', function(newCircleSizeType) {

              if (newCircleSizeType === undefined || scope.data === undefined) {
                  return;
              }

              var root  = scope.data.trees[scope.treeType],
                nodes = flatten(root);

              nodeDOM = vis.selectAll('circle').data(nodes, function(d) { return d.name; });

              var extent = d3.extent(nodes, function(d) { if (d.name[0] === 'B') { return d.size } })

              var size_lowerbound = 4;
              var sizeScale = d3.scale.linear()
                .domain(extent)
                .range([size_lowerbound, 10]);

              // transition from old to new
              nodeDOM.filter(function(d) { return d.name[0] === 'B' ? this : null; })
                .transition().delay(100).duration(750)
                .attr('r', function(d) { if (extent[0] == extent[1]) { return size_lowerbound; } else { return sizeScale(d.size);}});
          });

          scope.$watch('circleBorderType', function(newCircleBorderType) {

              if (newCircleBorderType === undefined || scope.data === undefined) {
                  return;
              }

              var root, nodes, colorExtent, colorScale, colorbarB;

              root  = scope.data.trees[scope.treeType];
              nodes = flatten(root);

              nodeDOM = vis.selectAll('circle').data(nodes, function(d) { return d.name; })

              if (scope.circleBorderType === 'SLogP') {
                borderScale = slogpScale;
                colorbarB = colorBar()
                  .color(borderScale).size(350).lineWidth(20).precision(4).tickFormat(d3.format('g'));
              } else if (scope.circleBorderType === 'Lig_Eff') {
                borderScale = ligeffScale;
                colorbarB = colorBar()
                  .color(borderScale).size(350).lineWidth(20).precision(4).tickFormat(d3.format('g'));
              } else {
                colorExtent = d3.extent(nodes, function(d) { if (d.name[0] === 'B') { return d.stroke; } });
                borderScale = d3.scale.linear()
                  .domain([colorExtent[0], d3.mean(colorExtent), colorExtent[1]])
                  .range(['#008000', '#FFFF00', '#FF0000']);
              }

              if (newCircleBorderType === 'None') {
                d3.select('.colorbarB').style('visibility','hidden');
              } else {
                d3.select('.colorbarB').style('visibility','visible').call(colorbarB);
              }

              // transition from old to new
              nodeDOM.filter(function(d) { return d.name[0] === 'B' ? this : null; })
                .transition().delay(100).duration(750)
                .style('stroke', function(d) { return borderScale(d.stroke); })
                .style('stroke-width', function(d) { return d.strokeWidth; });
          });

          scope.$watch('activityType', function(newActivityType) {

              if (newActivityType === undefined || scope.data === undefined) {
                  return;
              }

              var root  = scope.data.trees[scope.treeType],
                nodes = flatten(root);

              nodeDOM = vis.selectAll('circle').data(nodes, function(d) { return d.name; });

/*            var colorExtent = d3.extent(nodes, function(d) { if (d.name[0] === 'B') { return d.fill; } });
            activityScale = d3.scale.linear()
                .domain([colorExtent[0], d3.mean(colorExtent), colorExtent[1]])
                .range(['#008000', '#FFFF00', '#FF0000']);*/

              // transition from old to new
              nodeDOM.filter(function(d) { return d.name[0] === 'B' ? this : null; })
                .transition().delay(100).duration(750)
                .style('fill', function(d) { return activityScale(d.fill); });

          });

          scope.$watch('gravityValue', function(newGravity) {

              if (newGravity === undefined || scope.data === undefined) {
                  return;
              }

              //console.log('new gravity:' + newGravity);

              var force = scope.data.forces[scope.treeType];

              force.resume();
              force.gravity( newGravity );

          });

          scope.$watch('linkStrengthValue', function(newLinkStrength) {

              if (newLinkStrength === undefined || scope.data === undefined) {
                  return;
              }

              console.log('new linkStrength:' + newLinkStrength);

              var force = scope.data.forces[scope.treeType];

              force.linkStrength( newLinkStrength );
              force.resume();
              scope.$apply();

              console.log('tree type is: ', scope.treeType);
              console.log('setted linkStrength:' + force.linkStrength());
          });


          //watch the selected scope variable - this can be controlled from outside the directive or inside the directive
          scope.$watch('selected', function (selected) {

              if (selected === undefined || scope.data === undefined) {
                return;
              }

/*              var root  = scope.data.trees[scope.treeType],
                nodes = flatten(root);

              activityScale = d3.scale.linear()
                .domain(d3.extent(nodes, function(d) { return d.fill; }))
                .interpolate(d3.interpolateHcl)
                .range(['#008000', '#FFFF00', '#FF0000']);*/

              //if nothing is selected set the elements to not be selected
              if (selected === null) {
                nodeDOM
                  .filter(function(d) { return d.name[0] === 'B' ? this : null; })
                  .style('fill', function (d) {
                    return activityScale(d.fill);
                  });
              } else {
                nodeDOM
                  .filter(function(d) { return d.name[0] === 'B' ? this : null; })
                  .style('fill', function (d) {
                    return d.name === selected.id ? 'black' : activityScale(d.fill);
                  });
              }
            }
          );

      }
      return {
          restrict: 'E',
          link: link,
          scope: {
              data: '=',
              current: '=',
              gravityValue: '=',
              linkStrengthValue: '=',
              treeType: '=',
              circleSizeType: '=',
              circleBorderType: '=',
              activityType: '=',
              selected: '='
          }
      };
  });
