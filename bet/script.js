// Precomputed airfoil database: {airfoil: {Re: {alpha_key: [[alpha, Cl, Cd], ...]}}}
const airfoilDatabase = {
    "NACA0012": {
        "100000": {
            "-10 10 0.5": [
                [-10, -0.5, 0.02], [-9.5, -0.45, 0.019], [-9, -0.4, 0.018],
                [-8.5, -0.35, 0.017], [-8, -0.3, 0.016], [-7.5, -0.25, 0.015],
                [-7, -0.2, 0.014], [-6.5, -0.15, 0.013], [-6, -0.1, 0.012],
                [-5.5, -0.05, 0.011], [-5, 0.0, 0.01], [-4.5, 0.05, 0.01],
                [-4, 0.1, 0.01], [-3.5, 0.15, 0.01], [-3, 0.2, 0.01],
                [-2.5, 0.25, 0.01], [-2, 0.3, 0.01], [-1.5, 0.35, 0.01],
                [-1, 0.4, 0.01], [-0.5, 0.45, 0.01], [0, 0.5, 0.01],
                [0.5, 0.55, 0.01], [1, 0.6, 0.01], [1.5, 0.65, 0.01],
                [2, 0.7, 0.01], [2.5, 0.75, 0.01], [3, 0.8, 0.01],
                [3.5, 0.85, 0.01], [4, 0.9, 0.01], [4.5, 0.95, 0.01],
                [5, 1.0, 0.01], [5.5, 1.05, 0.011], [6, 1.1, 0.012],
                [6.5, 1.15, 0.013], [7, 1.2, 0.014], [7.5, 1.25, 0.015],
                [8, 1.3, 0.016], [8.5, 1.35, 0.017], [9, 1.4, 0.018],
                [9.5, 1.45, 0.019], [10, 1.5, 0.02]
            ]
        },
        "200000": {
            "-10 10 0.5": [
                [-10, -0.48, 0.018], [-9.5, -0.43, 0.017], [-9, -0.38, 0.016],
                [-8.5, -0.33, 0.015], [-8, -0.28, 0.014], [-7.5, -0.23, 0.013],
                [-7, -0.18, 0.012], [-6.5, -0.13, 0.011], [-6, -0.08, 0.01],
                [-5.5, -0.03, 0.009], [-5, 0.02, 0.008], [-4.5, 0.07, 0.008],
                [-4, 0.12, 0.008], [-3.5, 0.17, 0.008], [-3, 0.22, 0.008],
                [-2.5, 0.27, 0.008], [-2, 0.32, 0.008], [-1.5, 0.37, 0.008],
                [-1, 0.42, 0.008], [-0.5, 0.47, 0.008], [0, 0.52, 0.008],
                [0.5, 0.57, 0.008], [1, 0.62, 0.008], [1.5, 0.67, 0.008],
                [2, 0.72, 0.008], [2.5, 0.77, 0.008], [3, 0.82, 0.008],
                [3.5, 0.87, 0.008], [4, 0.92, 0.008], [4.5, 0.97, 0.008],
                [5, 1.02, 0.008], [5.5, 1.07, 0.009], [6, 1.12, 0.01],
                [6.5, 1.17, 0.011], [7, 1.22, 0.012], [7.5, 1.27, 0.013],
                [8, 1.32, 0.014], [8.5, 1.37, 0.015], [9, 1.42, 0.016],
                [9.5, 1.47, 0.017], [10, 1.52, 0.018]
            ]
        }
    },
    "NACA4412": {
        "100000": {
            "-10 10 0.5": [
                [-10, -0.6, 0.025], [-9.5, -0.55, 0.024], [-9, -0.5, 0.023],
                [-8.5, -0.45, 0.022], [-8, -0.4, 0.021], [-7.5, -0.35, 0.02],
                [-7, -0.3, 0.019], [-6.5, -0.25, 0.018], [-6, -0.2, 0.017],
                [-5.5, -0.15, 0.016], [-5, -0.1, 0.015], [-4.5, -0.05, 0.014],
                [-4, 0.0, 0.013], [-3.5, 0.05, 0.012], [-3, 0.1, 0.011],
                [-2.5, 0.15, 0.01], [-2, 0.2, 0.01], [-1.5, 0.25, 0.01],
                [-1, 0.3, 0.01], [-0.5, 0.35, 0.01], [0, 0.4, 0.01],
                [0.5, 0.45, 0.01], [1, 0.5, 0.01], [1.5, 0.55, 0.01],
                [2, 0.6, 0.01], [2.5, 0.65, 0.01], [3, 0.7, 0.01],
                [3.5, 0.75, 0.01], [4, 0.8, 0.01], [4.5, 0.85, 0.01],
                [5, 0.9, 0.01], [5.5, 0.95, 0.011], [6, 1.0, 0.012],
                [6.5, 1.05, 0.013], [7, 1.1, 0.014], [7.5, 1.15, 0.015],
                [8, 1.2, 0.016], [8.5, 1.25, 0.017], [9, 1.3, 0.018],
                [9.5, 1.35, 0.019], [10, 1.4, 0.02]
            ]
        }
    }
};

// Default airfoil data
const defaultAirfoilData = airfoilDatabase["NACA0012"]["100000"]["-10 10 0.5"];

// Linear interpolation function
function interp1(xTable, yTable, x) {
    const minX = Math.min(...xTable);
    const maxX = Math.max(...xTable);
    x = Math.max(minX, Math.min(maxX, x));

    for (let i = 0; i < xTable.length - 1; i++) {
        if (x >= xTable[i] && x <= xTable[i + 1]) {
            const x0 = xTable[i], x1 = xTable[i + 1];
            const y0 = yTable[i], y1 = yTable[i + 1];
            return y0 + (y1 - y0) * (x - x0) / (x1 - x0);
        }
    }
    return yTable[yTable.length - 1];
}

// Convert degrees to radians
function deg2rad(deg) {
    return deg * Math.PI / 180;
}

// Convert radians to degrees
function rad2deg(rad) {
    return rad * 180 / Math.PI;
}

// Generate linearly spaced array
function linspace(start, end, n) {
    const arr = [];
    const step = (end - start) / (n - 1);
    for (let i = 0; i < n; i++) {
        arr.push(start + i * step);
    }
    return arr;
}

// Input validation
function validateInputs(c, R, B, V_inf, omega, N, Re, alpha_range) {
    const errors = [];
    if (c <= 0) errors.push("Chord length must be positive.");
    if (R <= 0) errors.push("Blade radius must be positive.");
    if (B < 1) errors.push("Number of blades must be at least 1.");
    if (V_inf <= 0) errors.push("Free stream velocity must be positive.");
    if (N < 1) errors.push("Number of elements must be at least 1.");
    if (Re <= 0) errors.push("Reynolds number must be positive.");
    try {
        const [start, end, step] = alpha_range.trim().split(/\s+/).map(parseFloat);
        if (isNaN(start) || isNaN(end) || isNaN(step)) {
            errors.push("Alpha range must be three numbers (start end step).");
        } else if (step <= 0) {
            errors.push("Alpha step must be positive.");
        } else if (start >= end) {
            errors.push("Alpha start must be less than end.");
        }
    } catch {
        errors.push("Invalid alpha range format. Use: start end step (e.g., -10 10 0.5).");
    }
    return errors;
}

// Parse XFOIL file
async function parseXfoilFile(file) {
    const text = await file.text();
    const lines = text.split('\n').slice(12); // Skip 12 header lines, per MATLAB
    const data = lines
        .filter(line => line.trim())
        .map(line => {
            const [alpha, Cl, Cd] = line.trim().split(/\s+/).map(parseFloat);
            return [alpha, Cl, Cd];
        })
        .filter(row => !isNaN(row[0]) && !isNaN(row[1]) && !isNaN(row[2]));
    return data.length > 0 ? data : null;
}

document.getElementById('bet-form').addEventListener('submit', async function(event) {
    event.preventDefault();

    // Physical Parameters
    const rho = 1.225;
    const c = parseFloat(document.getElementById('chord').value);
    const R = parseFloat(document.getElementById('radius').value);
    const theta_inp = parseFloat(document.getElementById('twist').value);
    const theta = deg2rad(theta_inp);
    const B = parseInt(document.getElementById('blades').value);
    const V_inf = parseFloat(document.getElementById('velocity').value);
    const omega = parseFloat(document.getElementById('omega').value);
    const N = parseInt(document.getElementById('elements').value);
    const airfoil = document.getElementById('airfoil').value.trim().toUpperCase();
    const alpha_range = document.getElementById('alpha').value;
    const Re = parseInt(document.getElementById('reynolds').value);
    const xfoilFile = document.getElementById('xfoil-file').files[0];

    // Validate inputs
    const errors = validateInputs(c, R, B, V_inf, omega, N, Re, alpha_range);
    if (errors.length > 0) {
        document.getElementById('results').innerHTML = `
            <h3>Error</h3>
            <p>${errors.join('<br>')}</p>
        `;
        return;
    }

    // Parse alpha range
    const [alpha_start, alpha_end, alpha_step] = alpha_range.trim().split(/\s+/).map(parseFloat);
    const alpha_key = `${alpha_start} ${alpha_end} ${alpha_step}`;

    // Select airfoil data
    let selectedAirfoilData = defaultAirfoilData;
    let airfoilMessage = `Using default airfoil data (NACA0012, Re=100000, alpha=${alpha_key}).`;

    if (xfoilFile) {
        // Use uploaded XFOIL data
        const fileData = await parseXfoilFile(xfoilFile);
        if (fileData) {
            selectedAirfoilData = fileData;
            airfoilMessage = `Using uploaded XFOIL data from ${xfoilFile.name}.`;
        } else {
            airfoilMessage += `<br>Invalid XFOIL file. Using default data.`;
        }
    } else {
        // Use database
        const airfoilData = airfoilDatabase[airfoil];
        if (airfoilData) {
            const reData = airfoilData[Re.toString()];
            if (reData && reData[alpha_key]) {
                selectedAirfoilData = reData[alpha_key];
                airfoilMessage = `Using airfoil: ${airfoil}, Re=${Re}, alpha=${alpha_key}.`;
            } else {
                airfoilMessage = `No data for ${airfoil} at Re=${Re}, alpha=${alpha_key}. Using default.`;
            }
        } else {
            airfoilMessage = `Airfoil "${airfoil}" not found. Using default.`;
        }
    }

    // Load airfoil data
    const alpha_table = selectedAirfoilData.map(row => row[0]);
    const Cl_table = selectedAirfoilData.map(row => row[1]);
    const Cd_table = selectedAirfoilData.map(row => row[2]);

    // Discretization
    const r = linspace(0.1 * R, R, N);
    const dr = r[1] - r[0];

    let T = 0, Q = 0;
    const debugInfo = [];

    // Blade Element Theory calculations
    for (let i = 0; i < N; i++) {
        const V_a = V_inf;
        const V_t = omega * r[i];
        const V_res = Math.sqrt(V_a ** 2 + V_t ** 2);

        const phi = Math.atan2(V_a, V_t);
        const alpha = phi - theta;
        const alpha_deg = rad2deg(alpha);

        const Cl = interp1(alpha_table, Cl_table, alpha_deg);
        const Cd = interp1(alpha_table, Cd_table, alpha_deg);

        const dL = 0.5 * rho * (V_res ** 2) * c * Cl * dr;
        const dD = 0.5 * rho * (V_res ** 2) * c * Cd * dr;

        const dT = B * ((dL * Math.cos(phi)) - (dD * Math.sin(phi)));
        const dQ = B * r[i] * ((dL * Math.sin(phi)) + (dD * Math.cos(phi)));

        T += dT;
        Q += dQ;

        debugInfo.push(`Element ${i + 1}: r=${r[i].toFixed(2)}, V_res=${V_res.toFixed(2)}, alpha_deg=${alpha_deg.toFixed(2)}, Cl=${Cl.toFixed(2)}, Cd=${Cd.toFixed(2)}, dT=${dT.toFixed(2)}, dQ=${dQ.toFixed(2)}`);
    }

    const P = omega * Q;
    const eta = P !== 0 ? (T * V_inf) / P : 0;

    // Display results
    document.getElementById('results').innerHTML = `
        <h3>Simulation Results</h3>
        <p>${airfoilMessage}</p>
        <p>Total Thrust: ${T.toFixed(2)} N</p>
        <p>Total Torque: ${Q.toFixed(2)} Nm</p>
        <p>Power Required: ${P.toFixed(2)} W</p>
        <p>Prop Efficiency: ${eta.toFixed(2)}</p>
        <details>
            <summary>Debug Info</summary>
            <p style="font-size: 0.9em;">${debugInfo.join('<br>')}</p>
        </details>
    `;
});