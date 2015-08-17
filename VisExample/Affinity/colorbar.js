function colorBar(){
    var orient = "right",
    lineWidth = 40, //Function?... because that would be coooooool... not sure if it is compatible with axis.js
    size_ = 300,
    tickFormat = d3.format("3e"),
    color = d3.scale.linear().domain([0, 0.5, 1]).range(["blue", "green", "red"]), //v -> color
    line = d3.svg.line().interpolate("basis"),
    precision = 8,
    points_,
    tickSize_;

    function component(selection){
        selection.each(function(data,index){
            var container = d3.select(this),
            tickSize = tickSize_ || lineWidth,
            n,
            points = points_ || (((orient == "left") || (orient == "right"))?[[0,size_],[0,0]]:[[size_,0],[0,0]]),
            quads = quad(sample(line(points),precision)),
            size = (points)?n:size_,
            aScale = color.copy().interpolate(d3.interpolate).domain(d3.extent(color.domain())).range([size,0]), //v -> px
            colorExtent = d3.extent(color.domain()),
            normScale = color.copy().domain(color.domain().map(function(d){ return (d - colorExtent[0])/ (colorExtent[1] - colorExtent[0])})),

            //Save values for transitions
            oldLineWidth = this.__lineWidth__ || lineWidth;
            oldQuads = this.__quads__ || quads;
            this.__quads__ = quads;
            this.__lineWidth__ = lineWidth;

            //Enters
            var bar = container.selectAll("path.c").data(d3.range(quads.length), function(d){return d}),
            bEnter = bar.enter().insert("path","g.axis").classed("c",true),
            bExit = d3.transition(bar.exit()).remove(),
            bUpdate = d3.transition(bar),
            bTransform = function(selection,f,lw){
                selection.style("fill", function(d) { return normScale(f(d).t); })
                    .style("stroke", function(d) { return normScale(f(d).t); })
                    .attr("d", function(d) { var p = f(d); return lineJoin(p[0], p[1], p[2], p[3], lw); });};

            bEnter.call(bTransform,function(d){return oldQuads[oldQuads.length - 1]},oldLineWidth); // enter from last of oldQuad
            bExit.call(bTransform,function(d){return quads[quads.length - 1]},lineWidth); //exit from last of quads
            bUpdate.call(bTransform,function(d){return quads[d]},lineWidth)

            var colorBarAxis = d3.svg.axis().scale(aScale).orient(orient)
                .tickSize(tickSize).tickFormat(tickFormat),
            a = container.selectAll("g.axis").data(function(d){return (aScale)?[1]:[]}), //axis container
            aEnter = a.enter().append("g").classed("axis",true),
            aExit = d3.transition(a.exit()).remove(),
            aUpdate = d3.transition(a).call(colorBarAxis),
            aTransform = function(selection,lw){
                selection.attr("transform", "translate("
                               + (((orient == "right") || (orient == "left"))?-lw/2:0) + ","
                               + (((orient == "right") || (orient =="left"))?0:lw/2) + ")");};

            aEnter.call(aTransform,oldLineWidth);
            aExit.call(aTransform,lineWidth);
            aUpdate.call(aTransform,lineWidth);

            // Sample the SVG path string "d" uniformly with the specified precision.
            function sample(d,pre) {
                var path = document.createElementNS(d3.ns.prefix.svg, "path");
                path.setAttribute("d", d);

                n = path.getTotalLength();

                var t = [0], i = 0, dt = pre;
                while ((i += dt) < n) t.push(i);
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

            // Compute stroke outline for segment p12.
            function lineJoin(p0, p1, p2, p3, width) {
                var u12 = perp(p1, p2),
                r = width / 2,
                a = [p1[0] + u12[0] * r, p1[1] + u12[1] * r],
                b = [p2[0] + u12[0] * r, p2[1] + u12[1] * r],
                c = [p2[0] - u12[0] * r, p2[1] - u12[1] * r],
                d = [p1[0] - u12[0] * r, p1[1] - u12[1] * r];

                if (p0) { // clip ad and dc using average of u01 and u12
                    var u01 = perp(p0, p1), e = [p1[0] + u01[0] + u12[0], p1[1] + u01[1] + u12[1]];
                    a = lineIntersect(p1, e, a, b);
                    d = lineIntersect(p1, e, d, c);
                }

                if (p3) { // clip ab and dc using average of u12 and u23
                    var u23 = perp(p2, p3), e = [p2[0] + u23[0] + u12[0], p2[1] + u23[1] + u12[1]];
                    b = lineIntersect(p2, e, a, b);
                    c = lineIntersect(p2, e, d, c);
                }

                return "M" + a + "L" + b + " " + c + " " + d + "Z";
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

        })}

    component.orient = function(_) {
        if (!arguments.length) return orient;
        orient = _;
        return component;
    };

    component.lineWidth = function(_) {
        if (!arguments.length) return lineWidth;
        lineWidth = _;
        return component;
    };

    component.size = function(_) {
        if (!arguments.length) return size_;
        size_ = _;
        return component;
    };

    component.tickFormat = function(_) {
        if (!arguments.length) return tickFormat;
        tickFormat = _;
        return component;
    };

    component.tickSize = function(_) {
        if (!arguments.length) return tickSize_;
        tickSize_ = _;
        return component;
    };

    component.color = function(_) {
        if (!arguments.length) return color;
        color = _;
        return component;
    };

    component.precision = function(_) {
        if (!arguments.length) return precision;
        precision = _;
        return component;
    };

    component.points = function(_) {
        if (!arguments.length) return points_;
        points_ = _;
        return component;
    };

    component.line = function(_) {
        if (!arguments.length) return line;
        line = _;
        return component;
    };

    return component;
}
