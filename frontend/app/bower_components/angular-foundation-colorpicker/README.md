angular-foundation-colorpicker
==============================

Native AngularJS color picker directive for Foundation Zurb

<a href="http://web.hostdmk.net/github/colorpicker_foundation/" target="_blank">Demo page (Foundation 5)</a>

Installation
===============================
Copy css/colorpicker.css and js/foundation-colorpicker-module.js.
Add a dependency to your app, for instance:
angular.module('myApp', ['myApp.filters', 'myApp.services', 'myApp.directives', 'myApp.controllers', 'colorpicker.module'])

Examples:
===============================

Hex format
```html
<input colorpicker class="span2" type="text" ng-model="your_model" />
```
or
```html
<input colorpicker="hex" class="span2" type="text" ng-model="your_model" />
```

RGB format
```html
<input colorpicker="rgb" class="span2" type="text" ng-model="your_model" />
```

RBGA format
```html
<input colorpicker="rgba" class="span2" type="text" ng-model="your_model" />
```

As non input element
```html
<div colorpicker class="span2" ng-model="your_model"></div>
```

The color picker in a fixed element
```html
<input colorpicker colorpicker-fixed-position="true" class="span2" type="text" ng-model="your_model" />
```

When using fixed positioning, you can also put the picker into the parent element (this allows more styling control)
```html
<input colorpicker colorpicker-fixed-position="true" colorpicker-parent="true" class="span2" type="text" ng-model="your_model" />
```