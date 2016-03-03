/**
 * Created by ajing on 2/19/16.
 */

'use strict';

/**
 * @ngdoc service
 * @name frontendApp.settings
 * @description
 * # settings
 * Service in the frontendApp.
 */
angular.module('frontendApp')
  .service('settings', function () {

    this.defaultForce = {value: false};

    this.forceAct = {value: false};

  });