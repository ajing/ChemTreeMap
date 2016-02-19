'use strict';

/**
 * @ngdoc directive
 * @name frontendApp.directive:treeSlider
 * @description
 * # treeSlider
 */

angular.module('frontendApp')
  .directive('treeSlider', function () {
    return {
      scope: {
        min: '=',
        max: '=',
        value: '=',
        leftName: '=',
        rightName: '=',
        name: '=',
        id: '='
      },
      templateUrl: 'views/tree-slider.html',
      restrict: 'E',
      link: function postLink(scope, elements) {

        var nameTag = $('<div class="name section"></div>'),
          leftSection = $('<div class="left section"></div>'),
          leftLabel = $('<span class="sm-label">left</span>'),
          rightSection = $('<div class="right section"></div>'),
          rightLabel = $('<span class="sm-label">right</span>'),
          slider= $('<div class="slider"  title="Please activate the force directed graph first."></div>');

        //console.log('#' + scope.id);
        //console.log(slider);
        var el = elements.find('div');

        $(el).append([nameTag, slider]);

        var sliderScale = d3.scale.linear()
            .domain([scope.min, scope.max])
            .range([0, 280]);

        slider.append([leftSection, rightSection]);
        rightLabel.appendTo(rightSection);
        leftLabel.appendTo(leftSection);

        nameTag.text(scope.name);
        rightLabel.text(scope.rightName);
        leftLabel.text(scope.leftName);

        slider.slider({
          //range: true,
          min: scope.min,
          max: scope.max,
          step: 0.01,
          value: scope.value,
          //animate: 'fast',
          slide: function(evt, ui) {
              scope.value = ui.value;
              scope.$apply();
          }
        });

      }
    };
});
