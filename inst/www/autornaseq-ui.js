(function () {
  function ready(fn) {
    if (document.readyState === 'loading') {
      document.addEventListener('DOMContentLoaded', fn);
    } else {
      fn();
    }
  }

  function applyAdvancedClasses() {
    var panes = document.querySelectorAll('.content-wrapper .tab-pane');
    panes.forEach(function (pane) {
      pane.classList.add('arnx-pane');
      var boxes = pane.querySelectorAll('.box');
      boxes.forEach(function (box) {
        box.classList.add('arnx-pane-box');
      });
      if (boxes.length > 0) {
        boxes[0].classList.add('arnx-hero-box');
      }
    });
  }

  ready(function () {
    document.body.classList.add('autornaseq-loaded');
    applyAdvancedClasses();

    var menuLinks = document.querySelectorAll('.sidebar-menu a');
    menuLinks.forEach(function (link) {
      link.addEventListener('click', function () {
        document.body.classList.add('autornaseq-loaded');
        setTimeout(applyAdvancedClasses, 120);
      });
    });

    document.addEventListener('shown.bs.tab', function () {
      document.body.classList.add('autornaseq-loaded');
      setTimeout(applyAdvancedClasses, 120);
    });

    var target = document.querySelector('.content-wrapper');
    if (target && window.MutationObserver) {
      var observer = new MutationObserver(function () {
        applyAdvancedClasses();
      });
      observer.observe(target, { childList: true, subtree: true });
    }
  });
})();
