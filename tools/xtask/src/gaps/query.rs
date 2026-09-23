use super::model::{self, Capability, Issue, Ledger, Level, Result};

/// `list [--open] [--area <area>]`: one line per capability.
pub fn list(args: &[String]) -> Result<()> {
    let ledger = load()?;
    let open_only = args.iter().any(|a| a == "--open");
    let area = flag_value(args, "--area");
    for row in &ledger.capability {
        if open_only && row.level == Level::Implemented {
            continue;
        }
        if area.is_some_and(|area| row.area != area) {
            continue;
        }
        println!(
            "{:<4} {:<11} {:<15} {:<58} {}",
            row.id,
            row.level.as_str(),
            row.area,
            truncate(&row.name, 58),
            tracking(&ledger, row)
        );
    }
    Ok(())
}

/// `show <row-id | issue-key>`: everything needed to start work on it.
pub fn show(target: &str) -> Result<()> {
    let ledger = load()?;
    if let Some(row) = ledger.capability.iter().find(|row| row.id == target) {
        print_row(&ledger, row);
        return Ok(());
    }
    let issue = ledger.find_issue(target).ok_or_else(|| {
        format!("no row id, issue key or issue number `{target}`; try `cargo xtask gaps list`")
    })?;
    println!("{}  {}", issue.key, issue.title);
    println!("  tracked: {}", issue_ref(issue.number));
    println!("  priority: {}  effort: {}", issue.priority, issue.effort);
    for row in ledger
        .capability
        .iter()
        .filter(|row| row.issue.as_deref() == Some(issue.key.as_str()))
    {
        println!();
        print_row(&ledger, row);
    }
    Ok(())
}

/// Rows of `key` that are not yet implemented.
fn open_rows<'a>(ledger: &'a Ledger, key: &str) -> Vec<&'a str> {
    ledger
        .capability
        .iter()
        .filter(|row| row.issue.as_deref() == Some(key) && row.level != Level::Implemented)
        .map(|row| row.id.as_str())
        .collect()
}

/// `next`: open issues, highest priority first. An issue whose blocker
/// still has open rows is listed separately, so work starts where it can.
pub fn next() -> Result<()> {
    let ledger = load()?;
    let mut issues: Vec<_> = ledger.issue.iter().collect();
    issues.sort_by_key(|issue| {
        (
            rank(&issue.priority),
            rank(&issue.effort),
            issue.key.clone(),
        )
    });
    let mut blocked = Vec::new();
    println!("Ready now, highest priority and lowest effort first:");
    for issue in issues {
        let rows = open_rows(&ledger, &issue.key);
        if rows.is_empty() {
            continue;
        }
        let waiting: Vec<&str> = issue
            .blocked_by
            .iter()
            .filter(|key| !open_rows(&ledger, key).is_empty())
            .map(String::as_str)
            .collect();
        if !waiting.is_empty() {
            blocked.push((issue, waiting));
            continue;
        }
        print_issue_line(issue, &rows);
    }
    if !blocked.is_empty() {
        println!();
        println!("Blocked (start the blocker first):");
        for (issue, waiting) in blocked {
            println!(
                "  {:<27} {:<7} waits on {}",
                issue.key,
                issue_ref(issue.number),
                waiting.join(", ")
            );
        }
    }
    let untracked: Vec<&str> = ledger
        .capability
        .iter()
        .filter(|row| row.level != Level::Implemented && row.issue.is_none())
        .map(|row| row.id.as_str())
        .collect();
    println!();
    println!(
        "Open but untracked ({}): {}",
        untracked.len(),
        untracked.join(" ")
    );
    println!("  narrow rows that work for their named subset; file an issue before widening one.");
    println!();
    println!("Detail: cargo xtask gaps show <issue-key | row-id>");
    Ok(())
}

fn print_issue_line(issue: &Issue, rows: &[&str]) {
    println!(
        "  {:<27} {:<7} P={:<6} E={:<6} rows={}",
        issue.key,
        issue_ref(issue.number),
        issue.priority,
        issue.effort,
        rows.join(",")
    );
}

fn print_row(ledger: &Ledger, row: &Capability) {
    println!(
        "{} {} [{}, {}]",
        row.id,
        row.name,
        row.level.as_str(),
        row.area
    );
    println!("  {}", row.summary);
    println!("  tracked: {}", tracking(ledger, row));
    if !row.evidence.is_empty() {
        println!("  axiolid:");
        for cite in &row.evidence {
            println!("    {cite}");
        }
    }
    for (name, paths) in [("occt", &row.occt), ("cgal", &row.cgal)] {
        if paths.is_empty() {
            continue;
        }
        println!("  {name}:");
        if let Some(reference) = ledger.reference.get(name) {
            let tag = reference.tag.as_deref().unwrap_or("untagged");
            println!("    (pinned: {tag}, {})", &reference.commit[..12]);
        }
        for path in paths {
            println!("    {}", reference_url(ledger, name, path));
        }
    }
}

/// A link at the audited commit, so it keeps pointing at what was graded.
fn reference_url(ledger: &Ledger, library: &str, path: &str) -> String {
    match ledger.reference.get(library) {
        Some(reference) => format!(
            "{}/tree/{}/{path}",
            reference.repository.trim_end_matches('/'),
            reference.commit
        ),
        None => path.to_owned(),
    }
}

fn tracking(ledger: &Ledger, row: &Capability) -> String {
    match row.issue.as_deref().and_then(|key| ledger.issue(key)) {
        Some(issue) => format!("{} {}", issue_ref(issue.number), issue.key),
        None if row.level == Level::Implemented => "-".to_owned(),
        None => "untracked".to_owned(),
    }
}

fn issue_ref(number: u32) -> String {
    if number == 0 {
        "(unfiled)".to_owned()
    } else {
        format!("#{number}")
    }
}

fn rank(value: &str) -> u8 {
    match value {
        "Urgent" => 0,
        "High" => 1,
        "Medium" => 2,
        _ => 3,
    }
}

fn flag_value<'a>(args: &'a [String], flag: &str) -> Option<&'a str> {
    args.iter()
        .position(|a| a == flag)
        .and_then(|i| args.get(i + 1))
        .map(String::as_str)
}

fn truncate(text: &str, width: usize) -> String {
    if text.chars().count() <= width {
        text.to_owned()
    } else {
        text.chars().take(width - 1).chain(['…']).collect()
    }
}

fn load() -> Result<Ledger> {
    model::load(&model::workspace_root()?)
}
