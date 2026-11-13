(function () {
  Shiny.addCustomMessageHandler("cast_ui", function (state) {
    document.querySelectorAll("[data-cast-action]").forEach(function (button) {
      var available = true;
      if (button.id === "map_driver") available = state.map_available;
      if (button.id === "view" || button.id === "analyse") available = state.case_available;
      button.disabled = state.busy || !available;
    });
    var cancel = document.getElementById("cancel_job");
    if (cancel) cancel.hidden = !state.busy;
    ["bundle", "points", "figure"].forEach(function (name) {
      var link = document.getElementById("download_" + name);
      if (!link) return;
      var disabled = !state[name + "_available"];
      link.classList.toggle("disabled", disabled);
      link.setAttribute("aria-disabled", String(disabled));
      link.tabIndex = disabled ? -1 : 0;
    });
    var caseLink = document.getElementById("download_case");
    if (caseLink) {
      caseLink.classList.toggle("disabled", !state.case_export_available);
      caseLink.setAttribute("aria-disabled", String(!state.case_export_available));
      caseLink.tabIndex = state.case_export_available ? 0 : -1;
    }
  });
  Shiny.addCustomMessageHandler("cast_point", function (open) {
    var panel = document.getElementById("point_details");
    if (panel) panel.open = open;
  });
  document.addEventListener("click", function (event) {
    var button = event.target.closest("button[data-driver]");
    if (button) Shiny.setInputValue("driver_pick", button.dataset.driver, {priority: "event"});
  });
})();
