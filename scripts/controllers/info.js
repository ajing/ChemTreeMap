'use strict';

/**
 * @ngdoc function
 * @name frontendApp.controller:InfoCtrl
 * @description
 * # InfoCtrl
 * Controller of the frontendApp
 */
angular.module('frontendApp')
  .controller('InfoCtrl', function ($scope, $modalInstance) {

    $scope.dismiss = function () {
        $modalInstance.dismiss('cancel');
    };

});
