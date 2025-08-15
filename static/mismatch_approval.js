// mismatch_approval.js
//
// page loaded -> event listener for toggling mismatches & approving them; recalc percent and status
// 


// make sure everything runs once page is ready
document.addEventListener('DOMContentLoaded', function() {
  // add click 'event' to each mismatch (span)
  document.querySelectorAll('.mismatch').forEach(function(span) {
    span.addEventListener('click', function() {
       // toggle mismatch approval when click occurs
       span.classList.toggle('approved');
       recalculatePercentAndStatus(); // update both when toggle/click (inside event)
    });
  });
  recalculatePercentAndStatus();     // initial percent & status when page loaded
});

function recalculatePercentAndStatus() { // function for recalc
  const allSpans = document.querySelectorAll('.alignment-display span');
  let total = 0, matches = 0;
  // from total of alignments (match & mismatch); count mismatch as match if approved
  allSpans.forEach(span => {
    if (span.classList.contains('mismatch') || span.textContent === '|') {
        total++;
       if (span.classList.contains('approved') || span.textContent === '|') {
           matches++;
       }
    }
  });

  console.log('total:', total, 'matches:', matches); // de bug

  // calc percent identity
  const percent = total > 0 ? ((matches / total) * 100).toFixed(2) : 'N/A';
  document.getElementById('align_percent_display').textContent = percent + '%';

  // update hidden alignment_score for database_submission
  const alignScoreInput = document.getElementById('alignment_score');
  if (alignScoreInput) {
  alignScoreInput.value = percent;
  }
    
  // find arm and clone type (use re)
  function getArmAndClone(label) {
    const userInputDiv = document.querySelector('.user-input');
    if (!userInputDiv) return '';
    const regex = new RegExp(label + '\\s*(\\w+)', 'i');
    const match = userInputDiv.textContent.match(regex);
    return match ? match[1] : '';
  }
    
  // get arm and clone type; must be passed from results_display
  // it needs it to determine pass/fail (main.cgi has the pass/fail logic we draw from)
  const armType = getArmAndClone('Arm Type:');
  const cloneType = getArmAndClone('Clone Type:');
    
  let passStatus = 'Fail';
  if (cloneType === "Forward") {
    if (armType === "5Arm" && percent >= 95) passStatus = "Pass";
    if (armType === "3Arm" && percent >= 98) passStatus = "Pass";
  } else if (cloneType === "Reverse") {
    if (armType === "3Arm" && percent >= 95) passStatus = "Pass";
    if (armType === "5Arm" && percent >= 98) passStatus = "Pass";
  }
  // show pass status on page
  document.getElementById('pass_status_display').textContent = passStatus;
  
  //update hidden pass_status for database submission	
  const passStatusInput = document.getElementById('pass_status');
  if (passStatusInput) {
    passStatusInput.value = passStatus;
  }
}

// get another listener to submit form (even with live updates in beginning ^^, it must be recalc'd)
document.addEventListener('DOMContentLoaded', function() {
  var form = document.getElementById('SavedResultsForm');
  if (form) {
     // recalc once more
     form.addEventListener('submit', function() {
       recalculatePercentAndStatus(); 
     });
  }
});
