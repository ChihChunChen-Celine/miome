# ui.R
# -----------------------------------------------------------------------------
# Look and layout of the Shiny app: a navigation rail on the left, and for each
# analysis one page with a settings column next to the results. Settings and
# results scroll separately, so results stay in view while settings change.
# No web fonts or external files: the app works offline.
# -----------------------------------------------------------------------------

# One entry per page: title, colour of its zone, one-line description.
.PAGES <- list(
  data     = list(title = "Data",            col = "#8FB3D1",
                  desc = "Upload your files and check which samples were matched."),
  alpha    = list(title = "Alpha diversity", col = "#A9C59A",
                  desc = "Diversity within each sample, compared between groups."),
  pcoa     = list(title = "PCoA",            col = "#A9B4DE",
                  desc = "How similar the samples are to each other, as a map."),
  beta     = list(title = "Beta statistics", col = "#E79A9A",
                  desc = "Test whether groups differ in overall composition."),
  taxonomy = list(title = "Taxonomy",        col = "#E3BE6A",
                  desc = "Which taxa make up each sample, and a summary per group."),
  maaslin  = list(title = "MaAsLin3",        col = "#8FC3AC",
                  desc = "Find the individual taxa or pathways that differ between groups."),
  meas     = list(title = "Measurements",    col = "#D9A08B",
                  desc = "qPCR, metabolites and other measured values. Sequencing data is not needed.")
)

# Collapsible explanation, closed by default so the settings column stays short.
.help <- function(question, ...)
  shiny::tags$details(class = "mi-help", shiny::tags$summary(question), shiny::tags$div(...))

# Collapsible group of less-used settings.
.more <- function(label, ...)
  shiny::tags$details(class = "mi-more", shiny::tags$summary(label), shiny::tags$div(...))

.section <- function(label) shiny::tags$h3(class = "mi-section", label)
.step    <- function(n, label)
  shiny::tags$h3(class = "mi-section", shiny::tags$span(class = "mi-stepno", n), label)
.run_btn <- function(id, label) shiny::actionButton(id, label, class = "btn-primary btn-block mi-run")
.empty   <- function(...) shiny::tags$div(class = "mi-empty", ...)

# A white panel in the results area; `tools` (e.g. a download button) sit top right.
.card <- function(title = NULL, ..., tools = NULL, hint = NULL)
  shiny::tags$section(class = "mi-card",
    if (!is.null(title)) shiny::tags$div(class = "mi-card-head",
      shiny::tags$h2(title), if (!is.null(tools)) shiny::tags$div(class = "mi-tools", tools)),
    if (!is.null(hint)) shiny::tags$p(class = "mi-hint", hint),
    ...)

.page <- function(id, side, ..., active = FALSE) {
  z <- .PAGES[[id]]
  shiny::tags$div(id = paste0("page-", id), class = paste("mi-page", if (active) "on"),
    style = sprintf("--zone:%s", z$col),
    shiny::tags$header(class = "mi-head", shiny::tags$h1(z$title), shiny::tags$p(z$desc)),
    shiny::tags$div(class = "mi-side", side),
    shiny::tags$div(class = "mi-body", ...))
}

# Colour pickers, one per group, folded away until wanted. `cols` is a named
# vector: group -> the colour it has now.
.colour_editor <- function(cols, what = "group", max_n = 40) {
  if (!length(cols)) return(NULL)
  .more(sprintf("Choose a colour for each %s", what),
    if (length(cols) > max_n)
      shiny::tags$p(class = "mi-hint", sprintf("%d %ss is too many to colour one by one.", length(cols), what))
    else shiny::tagList(
      shiny::tags$div(class = "mi-cols", lapply(names(cols), function(l)
        shiny::tags$label(shiny::tags$input(type = "color", class = "mi-col", `data-level` = l,
                                            value = substr(cols[[l]], 1, 7)), l))),
      shiny::tags$p(class = "mi-hint", style = "margin:6px 0 0",
                    sprintf("A %s keeps its colour in every plot. ", what),
                    shiny::tags$a(href = "#", class = "mi-col-reset", "Back to the palette colours"))))
}
.grad_input <- function(id, label, value)
  shiny::tags$label(shiny::tags$input(type = "color", class = "mi-grad", id = id, value = value), label)

.nav_rail <- function() {
  link <- function(id, first = FALSE)
    shiny::tags$a(href = "#", `data-page` = id, class = if (first) "on",
                  `aria-current` = if (first) "page",
                  shiny::tags$span(class = "mi-dot", style = sprintf("background:%s", .PAGES[[id]]$col)),
                  .PAGES[[id]]$title)
  shiny::tags$aside(class = "mi-nav",
    shiny::tags$div(class = "mi-brand", "miome",
                    shiny::tags$small("microbiome analysis")),
    shiny::tags$nav(`aria-label` = "Analyses",
      link("data", TRUE),
      shiny::tags$div(class = "mi-navgroup", "Sequencing"),
      link("alpha"), link("pcoa"), link("beta"), link("taxonomy"), link("maaslin"),
      shiny::tags$div(class = "mi-navgroup", "Other data"),
      link("meas")),
    shiny::tags$div(class = "mi-navfoot",
      shiny::uiOutput("nav_status"),
      shiny::selectInput("palette", "Plot colours", choices = CAT_PALETTES, selected = "Miome"),
      shiny::checkboxInput("rev_pal", "Reverse order", FALSE)))
}

# Page switching happens in the browser. 'shown' / 'hidden' tell Shiny which
# outputs are visible, and the resize makes plots fit their new space.
.miome_js <- function() "
$(document).on('click', '.mi-nav a[data-page]', function(e) {
  e.preventDefault();
  var p = $(this).data('page');
  $('.mi-page').each(function() {
    var $e = $(this), on = this.id === 'page-' + p;
    if (on && !$e.hasClass('on')) { $e.trigger('show').addClass('on').trigger('shown'); }
    if (!on && $e.hasClass('on')) { $e.trigger('hide').removeClass('on').trigger('hidden'); }
  });
  $('.mi-nav a').removeClass('on').removeAttr('aria-current');
  $(this).addClass('on').attr('aria-current', 'page');
  $(window).trigger('resize');
});
// Colour pickers. One group keeps one colour in every plot, so a change is
// copied to the pickers of the same group on the other pages.
$(document).on('change', 'input.mi-col', function() {
  var l = $(this).attr('data-level'), v = this.value;
  $('input.mi-col').filter(function() { return $(this).attr('data-level') === l; }).val(v);
  Shiny.setInputValue('set_colour', {level: l, value: v}, {priority: 'event'});
});
$(document).on('change', 'input.mi-grad', function() { Shiny.setInputValue(this.id, this.value); });
$(document).on('click', 'a.mi-col-reset', function(e) {
  e.preventDefault(); Shiny.setInputValue('reset_colours', Math.random(), {priority: 'event'});
});"

.miome_css <- function() "
:root {
  --ground:#EEF1EA; --card:#FFFFFF; --ink:#2F3A36; --muted:#6D7873; --line:#DDE3DA;
  --rail:#F7F8F4; --action:#3F6F62; --action-dark:#33594F; --zone:#A9C59A;
  --serif:'Iowan Old Style','Palatino Linotype',Palatino,'Book Antiqua',Georgia,serif;
  --sans:'Avenir Next','Segoe UI Variable','Segoe UI',system-ui,-apple-system,sans-serif;
}
html, body { height:100%; }
body { background:var(--ground); color:var(--ink); font-family:var(--sans); font-size:13.5px;
       line-height:1.45; }
body > .container-fluid { padding:0; height:100%; }
.mi-app { display:flex; height:100vh; }

/* ---- navigation rail ---- */
.mi-nav { flex:0 0 204px; background:var(--rail); border-right:1px solid var(--line);
          display:flex; flex-direction:column; padding:22px 14px 14px; overflow-y:auto; }
.mi-brand { font-family:var(--serif); font-size:27px; line-height:1; padding:0 8px 20px;
            letter-spacing:.01em; }
.mi-brand small { display:block; font-family:var(--sans); font-size:11.5px; color:var(--muted);
                  margin-top:6px; letter-spacing:0; }
.mi-nav nav a { display:flex; align-items:center; gap:10px; padding:8px 10px; margin:1px 0;
                border-radius:9px; color:var(--ink); text-decoration:none; font-size:13.5px; }
.mi-nav nav a:hover { background:#ECEFE7; }
.mi-nav nav a.on { background:var(--card); box-shadow:0 0 0 1px var(--line); font-weight:600; }
.mi-nav nav a:focus-visible, .mi-help summary:focus-visible, .mi-more summary:focus-visible {
  outline:2px solid var(--action); outline-offset:1px; }
.mi-dot { width:11px; height:11px; border-radius:50%; flex:none; }
.mi-navgroup { font-size:11px; color:var(--muted); padding:14px 10px 4px; }
.mi-navfoot { margin-top:auto; padding:14px 8px 0; border-top:1px solid var(--line); }
.mi-navfoot .form-group { margin-bottom:4px; }
.mi-navfoot .checkbox { margin:2px 0 0; }
.mi-status { font-size:12px; color:var(--muted); margin-bottom:12px; }
.mi-status b { display:block; color:var(--ink); font-weight:600; }

/* ---- pages ---- */
.mi-main { flex:1; min-width:0; }
.mi-page { display:none; height:100vh; padding:0 22px 0 26px; column-gap:20px;
           grid-template-columns:340px minmax(0,1fr); grid-template-rows:auto minmax(0,1fr); }
.mi-page.on { display:grid; }
.mi-head { grid-column:1 / -1; padding:22px 0 14px; }
.mi-head h1 { font-family:var(--serif); font-size:28px; font-weight:500; margin:0; line-height:1.15; }
.mi-head h1::after { content:''; display:block; width:46px; height:4px; border-radius:2px;
                     background:var(--zone); margin-top:9px; }
.mi-head p { margin:9px 0 0; color:var(--muted); max-width:70ch; }
.mi-side, .mi-body { overflow-y:auto; padding-bottom:26px; }
.mi-side { background:var(--card); border-radius:14px 14px 0 0; border-top:4px solid var(--zone);
           padding:6px 18px 30px; }
.mi-body { padding-right:4px; }
.mi-card { background:var(--card); border-radius:14px; padding:16px 20px 18px; margin-bottom:16px; }
.mi-card-head { display:flex; align-items:center; justify-content:space-between; gap:12px;
                margin-bottom:8px; }
.mi-card h2 { font-family:var(--serif); font-size:18.5px; font-weight:500; margin:0; }
.mi-tools { display:flex; gap:8px; flex-wrap:wrap; }
.mi-hint { color:var(--muted); font-size:12.5px; margin:-2px 0 10px; max-width:80ch; }
.mi-section { font-family:var(--sans); font-size:13px; font-weight:700; margin:18px 0 10px;
              padding-top:14px; border-top:1px solid var(--line); display:flex; align-items:center;
              gap:8px; }
.mi-side > .mi-section:first-child, .mi-side > div:first-child > .mi-section:first-child {
  border-top:0; padding-top:4px; }
.mi-stepno { width:20px; height:20px; border-radius:50%; background:var(--zone); color:var(--ink);
             font-size:11.5px; display:inline-flex; align-items:center; justify-content:center; }
.mi-empty { color:var(--muted); padding:22px 4px; max-width:46ch; }
.mi-empty b { color:var(--ink); display:block; font-family:var(--serif); font-size:17px;
              font-weight:500; margin-bottom:4px; }

/* ---- form controls ---- */
label, .control-label { font-size:12.5px; font-weight:600; color:var(--ink); margin-bottom:4px; }
.checkbox label, .radio label, .radio-inline, .checkbox-inline { font-weight:400; font-size:13px; }
.form-group { margin-bottom:12px; }
.form-control, .selectize-input { border:1px solid #CDD5CB; border-radius:8px; box-shadow:none;
                                  font-size:13px; color:var(--ink); min-height:34px; }
.form-control:focus, .selectize-input.focus { border-color:var(--action);
                                              box-shadow:0 0 0 3px rgba(63,111,98,.16); }
.selectize-dropdown { border-radius:8px; border-color:#CDD5CB; font-size:13px; }
.selectize-dropdown .active { background:#E6EEE3; color:var(--ink); }
.selectize-control.multi .selectize-input > div { background:#E6EEE3; color:var(--ink);
                                                  border-radius:6px; padding:2px 6px; }
.btn { border-radius:8px; font-size:13px; }
.btn-default { background:var(--card); border-color:#CDD5CB; color:var(--ink); }
.btn-default:hover, .btn-default:focus { background:#F1F4EE; border-color:#B9C4B8; color:var(--ink); }
.btn-primary, .btn-primary:focus { background:var(--action); border-color:var(--action); color:#fff;
                                   font-weight:600; }
.btn-primary:hover, .btn-primary:active, .btn-primary:active:focus {
  background:var(--action-dark); border-color:var(--action-dark); }
.btn:focus-visible { outline:2px solid var(--action); outline-offset:2px; }
.mi-run { margin:14px 0 4px; padding:9px 12px; }
.btn-file { border-radius:8px 0 0 8px; }
.progress { border-radius:6px; height:6px; margin-top:4px; box-shadow:none; background:#E6EBE3; }
.progress-bar { background:var(--action); font-size:0; box-shadow:none; }
.shiny-input-container:not(.shiny-input-container-inline) { width:100%; max-width:none; }
hr { border-top-color:var(--line); }
pre { background:#F5F7F2; border:0; border-radius:10px; color:var(--ink); font-size:12.5px;
      padding:12px 14px; }
/* an empty result box must keep some size: Shiny does not fill outputs it sees as hidden */
pre.shiny-text-output:empty { visibility:hidden; height:1px; padding:0; margin:0; }

/* ---- notes ---- */
.info-box, .warn-box { border-radius:8px; padding:9px 12px; margin:8px 0; font-size:12.5px;
                       border-left:3px solid #8FB3D1; background:#EEF3F6; max-width:95ch; }
.warn-box { border-left-color:#D9A93F; background:#FBF2DD; }
.mi-help, .mi-more { margin:-4px 0 12px; font-size:12.5px; }
.mi-help summary, .mi-more summary { cursor:pointer; color:var(--action); border-radius:4px;
                                     width:fit-content; }
.mi-help > div { color:var(--muted); padding:6px 0 0 14px; max-width:60ch; }
.mi-more { margin:10px 0; }
.mi-more summary { font-weight:600; }
.mi-more > div { padding-top:10px; }

.mi-cols { display:flex; flex-wrap:wrap; gap:6px 16px; }
.mi-cols label { display:inline-flex; align-items:center; gap:7px; font-weight:400; font-size:13px;
                 margin:0; overflow-wrap:anywhere; }
input[type=color] { width:30px; height:24px; padding:0; border:1px solid #CDD5CB; border-radius:6px;
                    background:none; cursor:pointer; flex:none; }
.mi-col-reset { color:var(--action); }

/* ---- tables, tabs, messages ---- */
table.dataTable { font-size:12.5px; }
table.dataTable th, table.dataTable td { white-space:nowrap; }
table.dataTable thead th { border-bottom:1px solid var(--line) !important; font-weight:600; }
table.dataTable.no-footer { border-bottom:1px solid var(--line); }
.dataTables_wrapper .dataTables_paginate .paginate_button.current,
.dataTables_wrapper .dataTables_paginate .paginate_button.current:hover {
  background:#E6EEE3; border-color:transparent; border-radius:6px; }
.nav-pills > li > a { border-radius:8px; padding:6px 14px; color:var(--ink); }
.nav-pills > li.active > a, .nav-pills > li.active > a:hover, .nav-pills > li.active > a:focus {
  background:#E6EEE3; color:var(--ink); font-weight:600; }
.shiny-notification { border-radius:10px; border:1px solid var(--line); background:var(--card);
                      color:var(--ink); opacity:1; box-shadow:0 6px 24px rgba(47,58,54,.14); }
.shiny-output-error-validation { color:var(--muted); padding:18px 4px; }

@media (max-width: 1000px) {
  .mi-app { flex-direction:column; height:auto; }
  .mi-nav { flex:none; padding:12px; }
  .mi-nav nav { display:flex; flex-wrap:wrap; gap:2px; }
  .mi-navgroup { display:none; }
  .mi-brand { padding-bottom:10px; }
  .mi-navfoot { margin-top:10px; }
  .mi-page { height:auto; padding:0 14px; }
  .mi-page.on { display:block; }
  .mi-side { border-radius:14px; margin-bottom:16px; overflow:visible; }
  .mi-body { overflow:visible; }
}
@media (prefers-reduced-motion: no-preference) {
  .mi-nav nav a, .btn { transition:background-color .12s ease; }
}
"
