'use strict';


/**
 * @ngdoc overview
 * @name frontendApp
 * @description
 * # frontendApp
 *
 * Main module of the application.
 */

angular
  .module('frontendApp', [
    'oitozero.ngSweetAlert',
    'angular.filter',
    'colorpicker.module',
    'picardy.fontawesome',
    'mm.foundation',
    'angulartics',
    'angulartics.google.analytics',
    'ngAnimate',
    'ngCookies',
    'ngResource',
    'ngRoute',
    'ngSanitize',
    'ngTouch'
  ])
  .config(function ($routeProvider) {
    $routeProvider
      .when('/:dataset', {
          templateUrl: 'views/tree.html',
          controller: 'TreeController'
      })
      .otherwise({
        redirectTo: '/affinity'
      });
  });
