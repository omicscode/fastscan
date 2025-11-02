use anyhow::{Context, Result};
use bio::io::fasta;
use crossterm::{
    event::{self, DisableMouseCapture, EnableMouseCapture, Event, KeyCode},
    execute,
    terminal::{disable_raw_mode, enable_raw_mode, EnterAlternateScreen, LeaveAlternateScreen},
};
use ratatui::{
    backend::CrosstermBackend,
    layout::{Constraint, Direction, Layout},
    style::{Color, Modifier, Style},
    text::{Line, Span},
    widgets::{BarChart, Block, Borders, List, ListItem, ListState, Paragraph},
    Terminal,
};
use std::{
    env,
    fs::File,
    io::{self},
    path::{Path, PathBuf},
};
use walkdir::WalkDir;

#[derive(Clone)]
struct FastaFile {
    path: PathBuf,
    lengths: Vec<usize>,
}

struct App {
    files: Vec<FastaFile>,
    selected_file: usize,
    selected_seq: usize,
    lengths: Vec<usize>,
    bins: Vec<(String, u64)>,
    filter_min: usize,
    filter_max: usize,
    input: String,
    input_mode: bool,
    should_quit: bool,
}

impl Default for App {
    fn default() -> Self {
        Self {
            files: vec![],
            selected_file: 0,
            selected_seq: 0,
            lengths: vec![],
            bins: vec![],
            filter_min: 0,
            filter_max: usize::MAX,
            input: String::new(),
            input_mode: false,
            should_quit: false,
        }
    }
}

impl App {
    fn load_files_from_dir(&mut self, dir: &Path) -> Result<()> {
        self.files.clear();
        for entry in WalkDir::new(dir)
            .follow_links(true)
            .into_iter()
            .filter_map(|e| e.ok())
        {
            let path = entry.path();
            if path
                .extension()
                .and_then(|ext| ext.to_str())
                .map(|ext| ext == "fasta" || ext == "fa" || ext == "fna")
                .unwrap_or(false)
            {
                let reader = fasta::Reader::new(File::open(path)?);
                let mut lengths = Vec::new();
                for result in reader.records() {
                    let record = result?;
                    lengths.push(record.seq().len());
                }
                self.files.push(FastaFile {
                    path: path.to_path_buf(),
                    lengths,
                });
            }
        }
        self.files.sort_by(|a, b| a.path.cmp(&b.path));
        Ok(())
    }

    fn load_selected_file(&mut self) {
        if let Some(file) = self.files.get(self.selected_file) {
            self.lengths = file.lengths.clone();
            self.update_bins();
            self.apply_filter();
        }
    }

    fn update_bins(&mut self) {
        if self.lengths.is_empty() {
            self.bins = vec![];
            return;
        }
        let max_len = *self.lengths.iter().max().unwrap_or(&0);
        let bin_size = (max_len as f64 / 10.0).max(1.0) as usize;
        let num_bins = (max_len + bin_size - 1) / bin_size + 1;
        let mut bin_counts = vec![0u64; num_bins];

        for &len in &self.lengths {
            let bin_idx = (len / bin_size).min(num_bins - 1);
            bin_counts[bin_idx] += 1;
        }

        self.bins = bin_counts
            .into_iter()
            .enumerate()
            .map(|(i, count)| {
                (
                    format!("{}-{}", i * bin_size, (i + 1) * bin_size - 1),
                    count,
                )
            })
            .collect();
    }

    fn apply_filter(&mut self) {
        self.lengths
            .retain(|&len| len >= self.filter_min && len <= self.filter_max);
        self.update_bins();
    }

    fn next_file(&mut self) {
        if self.selected_file < self.files.len().saturating_sub(1) {
            self.selected_file += 1;
            self.load_selected_file();
        }
    }

    fn previous_file(&mut self) {
        if self.selected_file > 0 {
            self.selected_file -= 1;
            self.load_selected_file();
        }
    }

    fn next_seq(&mut self) {
        if self.selected_seq < self.lengths.len().saturating_sub(1) {
            self.selected_seq += 1;
        }
    }

    fn previous_seq(&mut self) {
        if self.selected_seq > 0 {
            self.selected_seq -= 1;
        }
    }

    fn set_filter_min(&mut self, min: usize) {
        self.filter_min = min;
        self.load_selected_file(); // Reapply from original
    }

    fn set_filter_max(&mut self, max: usize) {
        self.filter_max = max;
        self.load_selected_file();
    }
}

fn main() -> Result<()> {
    let args: Vec<String> = env::args().collect();
    let start_dir = if args.len() > 1 {
        Path::new(&args[1])
    } else {
        Path::new(".")
    };

    let mut app = App::default();
    app.load_files_from_dir(start_dir)
        .context("Failed to load FASTA files from directory")?;

    if app.files.is_empty() {
        eprintln!("No FASTA files found in directory.");
        return Ok(());
    }

    app.load_selected_file();

    // Setup terminal
    enable_raw_mode()?;
    let mut stdout = io::stdout();
    execute!(stdout, EnterAlternateScreen, EnableMouseCapture)?;
    let backend = CrosstermBackend::new(stdout);
    let mut terminal = Terminal::new(backend)?;

    let mut file_list_state = ListState::default();
    file_list_state.select(Some(app.selected_file));

    loop {
        terminal.draw(|f| ui(f, &mut app, &mut file_list_state))?;

        if let Event::Key(key) = event::read()? {
            if !app.input_mode {
                match key.code {
                    KeyCode::Char('q') => app.should_quit = true,
                    KeyCode::Down => app.next_file(),
                    KeyCode::Up => app.previous_file(),
                    KeyCode::Right => app.next_seq(),
                    KeyCode::Left => app.previous_seq(),
                    KeyCode::Char('m') => app.set_filter_min(1000),
                    KeyCode::Char('M') => app.set_filter_max(10000),
                    KeyCode::Char('i') => {
                        app.input_mode = true;
                        app.input.clear();
                    }
                    _ => {}
                }
            } else {
                // Input mode active
                match key.code {
                    KeyCode::Enter => {
                        let parts: Vec<&str> = app.input.trim().split('-').collect();
                        if parts.len() == 2 {
                            if let (Ok(min), Ok(max)) =
                                (parts[0].parse::<usize>(), parts[1].parse::<usize>())
                            {
                                app.set_filter_min(min);
                                app.set_filter_max(max);
                            }
                        }
                        app.input.clear();
                        app.input_mode = false;
                    }
                    KeyCode::Esc => {
                        app.input.clear();
                        app.input_mode = false;
                    }
                    KeyCode::Backspace if !app.input.is_empty() => {
                        app.input.pop();
                    }
                    KeyCode::Char(c) => {
                        app.input.push(c);
                    }
                    _ => {}
                }
            }

            if app.should_quit {
                break;
            }
        }
    }

    // Restore terminal
    disable_raw_mode()?;
    execute!(
        terminal.backend_mut(),
        LeaveAlternateScreen,
        DisableMouseCapture
    )?;
    terminal.show_cursor()?;

    Ok(())
}

fn ui(f: &mut ratatui::Frame, app: &mut App, file_list_state: &mut ListState) {
    let chunks = Layout::default()
        .direction(Direction::Vertical)
        .margin(2)
        .constraints([Constraint::Percentage(20), Constraint::Percentage(80)])
        .split(f.size());

    // === File List ===
    let file_items: Vec<ListItem> = app
        .files
        .iter()
        .enumerate()
        .map(|(i, file)| {
            let prefix = if i == app.selected_file { "> " } else { "  " };
            let style = if i == app.selected_file {
                Style::default()
                    .fg(Color::Yellow)
                    .add_modifier(Modifier::BOLD)
            } else {
                Style::default()
            };
            ListItem::new(format!("{}{}", prefix, file.path.display())).style(style)
        })
        .collect();

    let file_list = List::new(file_items)
        .block(Block::default().title("FASTA Files").borders(Borders::ALL))
        .highlight_style(Style::default().bg(Color::DarkGray))
        .highlight_symbol("> ");

    f.render_stateful_widget(file_list, chunks[0], file_list_state);

    // === Bottom Split ===
    let bottom_chunks = Layout::default()
        .direction(Direction::Horizontal)
        .constraints([Constraint::Percentage(50), Constraint::Percentage(50)])
        .split(chunks[1]);

    // === Histogram ===
    let histogram_widget: BarChart = if !app.bins.is_empty() {
        let max_count = app.bins.iter().map(|(_, c)| *c).max().unwrap_or(1);
        let bar_data: Vec<(&str, u64)> = app
            .bins
            .iter()
            .map(|(label, count)| (label.as_str(), count * 50 / max_count))
            .collect();

        BarChart::default()
            .data(&bar_data)
            .block(
                Block::default()
                    .title("Length Bins (filtered)")
                    .borders(Borders::ALL),
            )
            .bar_width(5)
            .bar_gap(1)
            .bar_style(Style::default().fg(Color::Green))
            .value_style(
                Style::default()
                    .fg(Color::Yellow)
                    .add_modifier(Modifier::BOLD),
            )
            .into()
    } else {
        BarChart::default()
    };
    f.render_widget(histogram_widget, bottom_chunks[0]);

    // === Right Panel: Filter + Info + Controls ===
    let filter_text = if app.input_mode {
        format!("Filter (min-max): {}", app.input)
    } else {
        format!(
            "Filter: {}-{} | Total: {} | Selected len: {}",
            app.filter_min,
            if app.filter_max == usize::MAX {
                "∞".to_string()
            } else {
                app.filter_max.to_string()
            },
            app.lengths.len(),
            app.lengths.get(app.selected_seq).copied().unwrap_or(0)
        )
    };

    let filter_para = Paragraph::new(filter_text)
        .block(
            Block::default()
                .title("Filter & Info")
                .borders(Borders::ALL),
        )
        .style(if app.input_mode {
            Style::default().fg(Color::Yellow)
        } else {
            Style::default()
        });

    let instructions = Paragraph::new(Line::from(vec![
        Span::styled("q", Style::default().fg(Color::Red)),
        Span::raw(":quit "),
        Span::styled("↑↓", Style::default().fg(Color::Cyan)),
        Span::raw(":files "),
        Span::styled("←→", Style::default().fg(Color::Cyan)),
        Span::raw(":seq "),
        Span::styled("i", Style::default().fg(Color::Magenta)),
        Span::raw(":filter "),
        Span::styled("m/M", Style::default().fg(Color::Green)),
        Span::raw(":1000/10000"),
    ]))
    .block(Block::default().title("Controls").borders(Borders::ALL));

    let right_chunks = Layout::default()
        .direction(Direction::Vertical)
        .constraints([Constraint::Percentage(60), Constraint::Percentage(40)])
        .split(bottom_chunks[1]);

    f.render_widget(filter_para, right_chunks[0]);
    f.render_widget(instructions, right_chunks[1]);
}
