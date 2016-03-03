/**
 * Created by ajing on 2/19/16.
 */
'use strict';

/**
 * @ngdoc function
 * @name frontendApp.controller:SettingsCtrl
 * @description
 * # SettingsCtrl
 * Controller of the frontendApp
 * Used to control the settings page
 */

angular.module('frontendApp')
  .controller('SettingsCtrl', function ($scope, $modalInstance, settings) {


    //dataService.current.forceAct = {value: false};
    $scope.forceAct = settings.forceAct;

    $scope.reset = function () {
      $scope.forceAct.value = settings.defaultForce.value;
    };

    $scope.dismiss = function () {
      $modalInstance.dismiss('cancel');
    };
  });