//
// function for coloring of target protein
//

function multiply(rgb1, rgb2) {
  var result = [],
  i = 0;
  for( ; i < rgb1.length; i++ ) {
    result.push(Math.floor(rgb1[i] * rgb2[i] / 255));
  }
  return result;
}
