'use strict';

/**
 * @ngdoc service
 * @name frontendApp.dataService
 * @description
 * # dataService
 * Factory in the frontendApp.
 */
angular.module('frontendApp')
  .factory('dataService', function ($http, $q) {
    // Service logic
    // ...

    // create the service to be returned by the factory
    var dataService = {};

    //selected molecules
    // also a variable dataService.model.selected
    dataService.model = {selected: null};

    //initially, there is no data
    dataService.data = null;

    //initially, there are no representations in data
    dataService.available = {
      treeTypes: [],
      circleSizeTypes: ['None'],
      circleBorderTypes: ['None'],
      activityTypes: []
    };

    //initially, there is no active spaces
    dataService.current = {
      treeType: null,
      circleSizeType: null,
      circleBorderType: null,
      activityType: null
    };

    dataService.flatten = function (root) {
      var nodes = [];

      function recurse(node) {
        if (node.children) { node.r = node.children.reduce(function(p, v) { return p + recurse(v); }, 0);}
        nodes.push(node);
        return node.size;
      }

      root.size = recurse(root);
      return nodes;
    };

    //
    // api for changing the data
    //

    // doesn't touch the data directly, calls setDimensionalityReductionType to actually change the data
    dataService.setTreeType = function (newTreeType) {

      console.log('Set Tree to ' + newTreeType);
      dataService.current.treeType = newTreeType;

      this.setActivityType(this.current.activityType);

    };

    //
    dataService.setCircleSizeType = function (newCircleSize) {

      console.log('Set circle size:' + newCircleSize);

      // should check if is a member of available
      this.current.circleSizeType = newCircleSize;

      var root  = this.data.trees[this.current.treeType];
      var nodes = dataService.flatten(root);

      if (newCircleSize === 'None') {
        nodes.forEach(function(d) {
          d.r = 2;
        });
      } else {
        nodes.forEach(function(d) {
          if (d.name.startsWith('B')) {
            d.r = dataService.data.compounds[parseInt(d.name.substring(1))].properties[newCircleSize];
          }
        });
      }

    };

    // use to change the border type
    dataService.setCircleBorderType = function (newCircleBorderType) {

      console.log('Set circleBorderType:', newCircleBorderType);

      // should check if is a member of available
      dataService.current.circleBorderType = newCircleBorderType;

      var root  = this.data.trees[this.current.treeType];
      var nodes = dataService.flatten(root);

      if (newCircleBorderType === 'None') {
        nodes.forEach(function(d) {
          d.stroke = 0;
          d.strokeWidth = 0;
        });
      } else {
        nodes.forEach(function(d) {
          if (d.name.startsWith('B')) {
            d.stroke = dataService.data.compounds[parseInt(d.name.substring(1))].properties[newCircleBorderType];
            d.strokeWidth = 3;
          }
        });
      }

    };

    // use to change the activity type
    dataService.setActivityType = function (newActivityType) {

      console.log('Set activity type:', newActivityType);

      // should check if is a member of available
      dataService.current.activityType = newActivityType;

      var root  = this.data.trees[this.current.treeType];
      var nodes = dataService.flatten(root);

      nodes.forEach(function(d) {
        if (d.name.startsWith('B')) {
          d.fill = dataService.data.compounds[parseInt(d.name.substring(1))].activities[newActivityType];
        }
      });

    };

    //remove all the data
    dataService.empty = function () {
      this.available.treeTypes.length = 0;
      this.available.circleSizeTypes = ['None'];
      this.available.circleBorderTypes = ['None'];
      this.available.activityTypes.length = 0;
      this.current.treeType = null;
      this.current.circleSizeType = null;
      this.current.circleBorderType = null;
      this.current.activityType = null;
      return this;
    };

    dataService.initializeData = function () {

      // remove old representation types if any previously existed
      this.empty();

      // add metadata for each option
      this.metadata = {};

      // add tree types
      this.data.metadata.treeTypes.forEach(function(d) {
        dataService.metadata[d.name] = d.metadata;
        dataService.available.treeTypes.push(d.name); });

      this.setTreeType(this.available.treeTypes[0]);

      // add circle size
      this.data.metadata.circleSizeTypes.forEach(function(d) {
        dataService.metadata[d.name] = d.metadata;
        dataService.available.circleSizeTypes.push(d.name); });

      this.setCircleSizeType(this.available.circleSizeTypes[0]);

      // add tree types
      this.data.metadata.circleBorderTypes.forEach(function(d) {
        dataService.metadata[d.name] = d.metadata;
        dataService.available.circleBorderTypes.push(d.name); });

      this.setCircleBorderType(this.available.circleBorderTypes[0]);

      // add activity types
      this.data.metadata.activityTypes.forEach(function(d) {
        dataService.metadata[d.name] = d.metadata;
        dataService.available.activityTypes.push(d.name); });

      this.setActivityType(this.available.activityTypes[0]);

      // set first in list to be active

      return this;
    };

    // load up an example dataset
    dataService.loadExample = function(name, callback) {

      var delay = $q.defer();

      $http.get('data/' + name + '.json')
        .then(function(response) {
          //retrieve the data as a property
          dataService.datasetName = name;
          dataService.data = response.data;

          //process the data
          dataService.initializeData();

          return delay.resolve(response);
        })
        .then(callback);

      return delay.promise;
    };

    return dataService;
  });
