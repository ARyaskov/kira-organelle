use std::path::Path;

use kira_organelle::run::derive_multi_report_dir_name;

#[test]
fn report_dir_name_from_input_dataset() {
    assert_eq!(
        derive_multi_report_dir_name(Path::new("data/GSE264586_RAW")),
        "GSE264586-report"
    );
    assert_eq!(
        derive_multi_report_dir_name(Path::new("data/GSE127465-RAW")),
        "GSE127465-report"
    );
    assert_eq!(
        derive_multi_report_dir_name(Path::new("data/experiment_set")),
        "experiment_set-report"
    );
}
