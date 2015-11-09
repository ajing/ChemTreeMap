/**
 * Created by ajing on 9/10/15.
 */


'use strict';

/**
 * @ngdoc function
 * @name frontendApp.controller:NavController
 * @description
 * # NavController
 * Controller of the frontendApp
 */
angular.module('frontendApp')
    .controller('NavController', function ($scope, $modal, dataService, SweetAlert) {

        $scope.dataService = dataService;

        $scope.currentSearch = '';


        $scope.getInfo = function(e, infoObject) {
            e.stopPropagation();
            //console.log(dataService.metadata);
            SweetAlert.swal({
                title: infoObject,
                html: dataService.metadata[infoObject],
                allowOutsideClick: true
            });
        };

        $scope.select = function(a) {
            dataService.model.selected = a;
            console.log(dataService.model.selected);
        };

        //open the info modal
        $scope.openInfo = function () {

            $modal.open({
                    templateUrl: 'views/info.html',
                    controller: 'InfoCtrl'
                }
            );
        };
    });
