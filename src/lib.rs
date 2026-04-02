pub mod consts;
pub mod mixtures;
pub mod properties;

#[macro_export]
macro_rules! assert_delta {
    ($left:expr, $right:expr, $delta:expr) => {{
        let diff = ($left - $right).abs();
        assert!(
            diff <= $delta,
            "|{} - {}| = {} > {}",
            $left,
            $right,
            diff,
            $delta
        );
    }};
}

#[macro_export]
macro_rules! assert_vec_delta {
    ($left:expr, $right:expr, $delta:expr) => {{
        let left = &$left;
        let right = &$right;
        assert_eq!(
            left.len(),
            right.len(),
            "slice lengths differ: {} != {}",
            left.len(),
            right.len()
        );
        for (i, (l, r)) in left.iter().zip(right.iter()).enumerate() {
            let diff = (*l - *r).abs();
            assert!(
                diff <= $delta,
                "element {} not within delta: |{} - {}| = {} > {}",
                i,
                l,
                r,
                diff,
                $delta
            );
        }
    }};
}
