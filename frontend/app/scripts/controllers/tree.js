'use strict';

/**
 * @ngdoc function
 * @name frontendApp.controller:TreeController
 * @description
 * # TreeController
 * Controller of the frontendApp
 */
angular.module('frontendApp')
  .controller('TreeController', function ($scope, $routeParams, dataService, settings) {

    //locate database from the route
    $scope.datasetName = $routeParams.dataset;

    //selected information is currently nothing
    $scope.model = dataService.model;

    $scope.current = dataService.current;

    $scope.flatten = dataService.flatten;

    $scope.tooltip = {visibility: false};

    $scope.settings = settings;

    $scope.$watch('tooltip.visibility', function(newVis) {
      if ( newVis === false ) {
        $scope.model.selected = null;
      }
    });

    $scope.select = function(d) {
      $scope.model.selected = d;
    };

    $scope.$watch('model.selected', function(selected) {
      console.log('selected:', selected);

      if (selected === null) {
        $scope.tooltip.visibility = false;
      } else {
        //set up the tooltip for the specific selected item
        $scope.tooltip.visibility = true;
        $scope.tooltip.data = {
          compound: true,
          object: selected
        };
      }
    });



    $scope.gravitySlider      = {min: 0, max: 0.2, value:0.1, id: 'gravity', name: 'Radius of Display', leftName: 'Lager', rightName: 'Smaller'};

    $scope.linkStrengthSlider = {min: 0, max: 10, value: 1, id: 'linkstrength', name: 'Compactness', leftName: 'Looser', rightName: 'Tight'};

    dataService.loadExample($routeParams.dataset, function() {
        $scope.data = dataService.data;
/*        // i deleted all slider parameters, but I need to add color bar parameters here.
        $scope.$watch('dataService.current.circleBorderType', function() {
            var extent = d3.extent($scope.data.nodes, function(d) { return d.strock; });
            $scope.borderColorbar.min = extent[0];
            $scope.borderColorbar.max = extent[1];
        });*/

    });
});
