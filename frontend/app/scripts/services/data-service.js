'use strict';

/**
 * @ngdoc service
 * @name frontendApp.dataService
 * @description
 * # dataService
 * Factory in the frontendApp.
 */
angular.module('frontendApp')
  .factory('dataService', function ($http, $q) {
    // Service logic
    // ...

    // create the service to be returned by the factory
    var dataService = {};

    //selected molecules
    // also a variable dataService.model.selected
    dataService.model = {selected: null};

    //initially, there is no data
    dataService.data = null;

    //initially, there are no representations in data
    dataService.available = {
      treeTypes: [],
      circleSizeTypes: ['None'],
      circleBorderTypes: ['None'],
      activityTypes: []
    };

    //initially, there is no active spaces
    dataService.current = {
      treeType: null,
      circleSizeType: null,
      circleBorderType: null,
      activityType: null
    };

    dataService.flatten = function (root) {
      var nodes = [];

      function recurse(node) {
        if (node.children) { node.r = node.children.reduce(function(p, v) { return p + recurse(v); }, 0);}
        nodes.push(node);
        return node.size;
      }

      root.size = recurse(root);
      return nodes;
    };

    //
    // api for changing the data
    //

    // doesn't touch the data directly, calls setDimensionalityReductionType to actually change the data
    dataService.setTreeType = function (newTreeType) {

      console.log('Set Tree to ' + newTreeType);
      dataService.current.treeType = newTreeType;

      this.setActivityType(this.current.activityType);

    };

    //
    dataService.setCircleSizeType = function (newCircleSize) {

      console.log('Set circle size:' + newCircleSize);

      // should check if is a member of available
      this.current.circleSizeType = newCircleSize;

      var root  = this.data.trees[this.current.treeType];
      var nodes = dataService.flatten(root);

      if (newCircleSize === 'None') {
        nodes.forEach(function(d) {
          d.r = 2;
        });
      } else {
        nodes.forEach(function(d) {
          if (d.name.startsWith('B')) {
            d.r = dataService.data.compounds[parseInt(d.name.substring(1))].properties[newCircleSize];
          }
        });
      }

    };

    // use to change the border type
    dataService.setCircleBorderType = function (newCircleBorderType) {

      console.log('Set circleBorderType:', newCircleBorderType);

      // should check if is a member of available
      dataService.current.circleBorderType = newCircleBorderType;

      var root  = this.data.trees[this.current.treeType];
      var nodes = dataService.flatten(root);

      if (newCircleBorderType === 'None') {
        nodes.forEach(function(d) {
          d.stroke = 0;
          d.strokeWidth = 0;
        });
      } else {
        nodes.forEach(function(d) {
          if (d.name.startsWith('B')) {
            d.stroke = dataService.data.compounds[parseInt(d.name.substring(1))].properties[newCircleBorderType];
            d.strokeWidth = 3;
          }
        });
      }

    };

    // use to change the activity type
    dataService.setActivityType = function (newActivityType) {

      console.log('Set activity type:', newActivityType);

      // should check if is a member of available
      dataService.current.activityType = newActivityType;

      var root  = this.data.trees[this.current.treeType];
      var nodes = dataService.flatten(root);

      nodes.forEach(function(d) {
        if (d.name.startsWith('B')) {
          d.fill = dataService.data.compounds[parseInt(d.name.substring(1))].activities[newActivityType];
        }
      });

    };

    //remove all the data
    dataService.empty = function () {
      this.available.treeTypes.length = 0;
      this.available.circleSizeTypes = ['None'];
      this.available.circleBorderTypes = ['None'];
      this.available.activityTypes.length = 0;
      this.current.treeType = null;
      this.current.circleSizeType = null;
      this.current.circleBorderType = null;
      this.current.activityType = null;
      return this;
    };

    dataService.initializeData = function () {

      // remove old representation types if any previously existed
      this.empty();

      // add metadata for each option
      this.metadata = {};

      // add tree types
      this.data.metadata.treeTypes.forEach(function(d) {
        dataService.metadata[d.name] = d.metadata;
        dataService.available.treeTypes.push(d.name); });

      this.setTreeType(this.available.treeTypes[0]);

      // add circle size
      this.data.metadata.circleSizeTypes.forEach(function(d) {
        dataService.metadata[d.name] = d.metadata;
        dataService.available.circleSizeTypes.push(d.name); });

      this.setCircleSizeType(this.available.circleSizeTypes[0]);

      // add tree types
      this.data.metadata.circleBorderTypes.forEach(function(d) {
        dataService.metadata[d.name] = d.metadata;
        dataService.available.circleBorderTypes.push(d.name); });

      this.setCircleBorderType(this.available.circleBorderTypes[0]);

      // add activity types
      this.data.metadata.activityTypes.forEach(function(d) {
        dataService.metadata[d.name] = d.metadata;
        dataService.available.activityTypes.push(d.name); });

      this.setActivityType(this.available.activityTypes[0]);

      // set first in list to be active

      return this;
    };

    // load up an example dataset
    dataService.loadExample = function(name, callback) {

      var delay = $q.defer();

      $http.get('data/' + name + '.json')
        .then(function(response) {
          //retrieve the data as a property
          dataService.datasetName = name;
          dataService.data = response.data;

          //process the data
          dataService.initializeData();

          return delay.resolve(response);
        })
        .then(callback);

      return delay.promise;
    };

    dataService.colorBar = function(){
      var orient = 'right',
        lineWidth = 40, //Function?... because that would be coooooool... not sure if it is compatible with axis.js
        size_ = 300,
        tickFormat = d3.format('3e'),
        color = d3.scale.linear().domain([0, 0.5, 1]).range(['blue', 'green', 'red']), //v -> color
        line = d3.svg.line().interpolate('basis'),
        precision = 8,
        points_,
        tickSize_,
        oldLineWidth,
        oldQuads;


      // Sample the SVG path string 'd' uniformly with the specified precision.
      function sample(d,pre) {
        var path = document.createElementNS(d3.ns.prefix.svg, 'path');
        path.setAttribute('d', d);

        var n = path.getTotalLength();

        var t = [0], i = 0, dt = pre;
        while ((i += dt) < n) {
          t.push(i);
        }
        t.push(n);

        return t.map(function(t) {
          var p = path.getPointAtLength(t), a = [p.x, p.y];
          a.t = t / n;
          return a;
        });

      }

      // Compute quads of adjacent points [p0, p1, p2, p3].
      function quad(pts) {
        return d3.range(pts.length - 1).map(function(i) {
          var a = [pts[i - 1], pts[i], pts[i + 1], pts[i + 2]];
          a.t = (pts[i].t + pts[i + 1].t) / 2;
          return a;
        });
      }

      function component(selection){
        selection.each(function(){
          var container = d3.select(this),
            tickSize = tickSize_ || lineWidth,
            n,
            points = points_ || (((orient === 'left') || (orient === 'right'))?[[0,size_],[0,0]]:[[size_,0],[10,0]]),
            quads = quad(sample(line(points),precision)),
            size = (points)?n:size_,
            aScale = color.copy().interpolate(d3.interpolate).domain(color.domain()).range([size,0]), //v -> px
            colorExtent = color.domain(),
            normScale = color.copy().domain(color.domain().map(function(d){ return (d - colorExtent[0])/ (colorExtent[1] - colorExtent[0]);}));

          //Save values for transitions
          oldLineWidth = this.__lineWidth__ || lineWidth;
          oldQuads = this.__quads__ || quads;
          this.__quads__ = quads;
          this.__lineWidth__ = lineWidth;

          //Enters
          var bar = container.selectAll('path.c').data(d3.range(quads.length), function(d){return d;}),
            bEnter = bar.enter().insert('path','g.axis').classed('c',true),
            bExit = d3.transition(bar.exit()).remove(),
            bUpdate = d3.transition(bar),
            bTransform = function(selection,f,lw){
              selection.style('fill', function(d) { return normScale(f(d).t); })
                .style('stroke', function(d) { return normScale(f(d).t); })
                .attr('d', function(d) { var p = f(d); return lineJoin(p[0], p[1], p[2], p[3], lw); });};

          bEnter.call(bTransform,function(){return oldQuads[oldQuads.length - 1];},oldLineWidth); // enter from last of oldQuad
          bExit.call(bTransform,function(){return quads[quads.length - 1];},lineWidth); //exit from last of quads
          bUpdate.call(bTransform,function(d){return quads[d];},lineWidth);

          var colorBarAxis = d3.svg.axis().scale(aScale).orient(orient)
              .tickSize(tickSize).tickFormat(tickFormat),
            a = container.selectAll('g.axis').data(function(){return (aScale)?[1]:[];}), //axis container
            aEnter = a.enter().append('g').classed('axis',true),
            aExit = d3.transition(a.exit()).remove(),
            aUpdate = d3.transition(a).call(colorBarAxis),
            aTransform = function(selection,lw){
              selection.attr('transform', 'translate(' + (((orient === 'right') || (orient === 'left'))?-lw/2:0) + ',' + (((orient === 'right') || (orient ==='left'))?0:lw/2) + ')');};

          aEnter.call(aTransform,oldLineWidth);
          aExit.call(aTransform,lineWidth);
          aUpdate.call(aTransform,lineWidth);

          // Compute stroke outline for segment p12.
          function lineJoin(p0, p1, p2, p3, width) {
            var u12 = perp(p1, p2), e, u23,
              r = width / 2,
              a = [p1[0] + u12[0] * r, p1[1] + u12[1] * r],
              b = [p2[0] + u12[0] * r, p2[1] + u12[1] * r],
              c = [p2[0] - u12[0] * r, p2[1] - u12[1] * r],
              d = [p1[0] - u12[0] * r, p1[1] - u12[1] * r];

            if (p0) { // clip ad and dc using average of u01 and u12
              var u01 = perp(p0, p1);
              e = [p1[0] + u01[0] + u12[0], p1[1] + u01[1] + u12[1]];
              a = lineIntersect(p1, e, a, b);
              d = lineIntersect(p1, e, d, c);
            }

            if (p3) { // clip ab and dc using average of u12 and u23
              u23 = perp(p2, p3);
              e = [p2[0] + u23[0] + u12[0], p2[1] + u23[1] + u12[1]];
              b = lineIntersect(p2, e, a, b);
              c = lineIntersect(p2, e, d, c);
            }

            return 'M' + a + 'L' + b + ' ' + c + ' ' + d + 'Z';
          }

          // Compute intersection of two infinite lines ab and cd.
          function lineIntersect(a, b, c, d) {
            var x1 = c[0], x3 = a[0], x21 = d[0] - x1, x43 = b[0] - x3,
              y1 = c[1], y3 = a[1], y21 = d[1] - y1, y43 = b[1] - y3,
              ua = (x43 * (y1 - y3) - y43 * (x1 - x3)) / (y43 * x21 - x43 * y21);
            return [x1 + ua * x21, y1 + ua * y21];
          }

          // Compute unit vector perpendicular to p01.
          function perp(p0, p1) {
            var u01x = p0[1] - p1[1], u01y = p1[0] - p0[0],
              u01d = Math.sqrt(u01x * u01x + u01y * u01y);
            return [u01x / u01d, u01y / u01d];
          }

        });}

      component.orient = function(_) {
        if (!arguments.length) {
          return orient;
        }
        orient = _;
        return component;
      };

      component.lineWidth = function(_) {
        if (!arguments.length) {
          return lineWidth;
        }
        lineWidth = _;
        return component;
      };

      component.size = function(_) {
        if (!arguments.length) {
          return size_;
        }
        size_ = _;
        return component;
      };

      component.tickFormat = function(_) {
        if (!arguments.length) {
          return tickFormat;
        }
        tickFormat = _;
        return component;
      };

      component.tickSize = function(_) {
        if (!arguments.length) {
          return tickSize_;
        }
        tickSize_ = _;
        return component;
      };

      component.color = function(_) {
        if (!arguments.length) {
          return color;
        }
        color = _;
        return component;
      };

      component.precision = function(_) {
        if (!arguments.length) {
          return precision;
        }
        precision = _;
        return component;
      };

      component.points = function(_) {
        if (!arguments.length) {
          return points_;
        }
        points_ = _;
        return component;
      };

      component.line = function(_) {
        if (!arguments.length) {
          return line;
        }
        line = _;
        return component;
      };

      return component;
    }


    return dataService;
  });
