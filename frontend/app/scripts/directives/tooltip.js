'use strict';

/**
 * @ngdoc directive
 * @name frontendApp.directive:tooltip
 * @description
 * # tooltip
 */
angular.module('frontendApp')
  .directive('tooltips', function (dataService) {

    return {
        scope: {visibility: '=', data: '=', height: '@', click: '='},
        templateUrl: 'views/tooltip.html',
        restrict: 'E',
        link: function (scope, elements) {

            var parent = $(elements[0]);

            scope.dataService = dataService;

            scope.close = function() {
                scope.visibility = false;
            };


            scope.$watch('visibility', function(newVisibility) {
                if (newVisibility) {
                    parent.show(400);
                } else {
                    parent.hide(400);
                }
            });


        }
    };
  });
