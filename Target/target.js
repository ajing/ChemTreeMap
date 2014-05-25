var tproteins = ["first", "second", "third"];

var targetprotein = (function(group){
  // create table
  var table = d3.select('.stat').selectAll('table')
    .data([group])
  .enter()
    .append('table')
    .attr('tabindex', 1); // enable focus

  var tr = table.selectAll('tr')
    .data(function(d) { return d; })
  .enter()
    .append('tr')
    .classed('selected', function(d) { return d._selected; })

  tr.append('td').text(function(d, i) { return i + 1; });
  tr.append('td').text(function(d) { return d; });


  // color mixing for properties
  function multiply(rgb1, rgb2) {
    var result = [],
    i = 0;
    for( ; i < rgb1.length; i++ ) {
      result.push(Math.floor(rgb1[i] * rgb2[i] / 255));
    }
    return result;
  }

  // Update color
  var colortrans = d3.scale.ordinal().domain(tproteins)
      .range(colorbrewer.Spectral[tproteins.length])
  function updatecolor(selected){
      node.style("fill", function(d){
        if ( tproteins.indexOf(d.target) != -1) {
          return colortrans(d.target);
        }
        return d.group;
      })
  }

  // JQuery
  $(function() {
    $( ".stat" ).selectable({
      filter:'tr',
      stop: function() {
        var result = [];
        $( ".ui-selected", this ).each(function() {
          var index = $( ".stat tr" ).index( this );
          result.push( group[index] );
        });
        console.log(result);
        updatecolor(result);
      }
    });
  });
}(tproteins))
