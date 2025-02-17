// Course data
const courses = [
    {"id": "חדוא 1", "name": {"he": "חדוא 1", "en": "Calculus 1"}, "prereqs": [], "coreqs": [], "category": "analysis"},
    {"id": "לינארית 1", "name": {"he": "לינארית 1", "en": "Linear Algebra 1"}, "prereqs": [], "coreqs": [], "category": "algebra"},
    {"id": "מבוא לתורת הקבוצות", "name": {"he": "מבוא לתורת הקבוצות", "en": "Introduction to Set Theory"}, "prereqs": [], "coreqs": [], "category": "foundations"},
    {"id": "מבוא מורחב למדמ״ח", "name": {"he": "מבוא מורחב למדמ״ח", "en": "Extended Introduction to Computer Science"}, "prereqs": [], "coreqs": ["מבוא לתורת הקבוצות"], "category": "computer_science"},
    {"id": "קומבינטוריקה", "prereqs": ["חדוא 1", "לינארית 1"], "coreqs": ["מבוא לתורת הקבוצות"], "category": "discrete"},
    {"id": "חדוא 2", "prereqs": ["חדוא 1"], "coreqs": [], "category": "analysis"},
    {"id": "לינארית 2", "prereqs": ["לינארית 1"], "coreqs": [], "category": "algebra"},
    {"id": "מבוא להסתברות", "prereqs": ["קומבינטוריקה","מבוא לתורת הקבוצות"], "coreqs": ["חדוא 2","לינארית 2"], "category": "probability"},
    {"id": "חבורות", "prereqs": [], "coreqs": ["לינארית 2"], "category": "algebra"},
    {"id": "קומבינטוריקה בסיסית", "prereqs": ["קומבינטוריקה"], "coreqs": [], "category": "discrete"},
    {"id": "תורת הגרפים", "prereqs": ["מבוא להסתברות"], "coreqs": [], "category": "discrete"},
    {"id": "שיטות מתמטיות בתורת המשחקים", "prereqs": ["מבוא להסתברות"], "coreqs": [], "category": "applications"},
    {"id": "תורת המספרים", "prereqs": ["חדוא 2", "לינארית 2"], "coreqs": [], "category": "algebra"},
    {"id": "חדוא 3", "prereqs": ["חדוא 2", "לינארית 2", "קומבינטוריקה"], "coreqs": [], "category": "analysis"},
    {"id": "תורת ההצגות", "prereqs": ["חבורות"], "coreqs": [], "category": "algebra"},
    {"id": "גלואה", "prereqs": ["חבורות"], "coreqs": [], "category": "algebra"},
    {"id": "אלגברה קומוטטיבית", "prereqs": ["חבורות"], "coreqs": [], "category": "algebra"},
    {"id": "ODE", "prereqs": ["לינארית 2"], "coreqs": ["חדוא 3"], "category": "analysis"},
    {"id": "גאומטריה דיפרנציאלית", "prereqs": ["לינארית 2"], "coreqs": ["חדוא 3", "ODE"], "category": "geometry"},
    {"id": "ממשיות", "prereqs": ["לינארית 2"], "coreqs": ["חדוא 3"], "category": "analysis"},
    {"id": "מרוכבות 1", "prereqs": ["חדוא 2"], "coreqs": ["חדוא 3"], "category": "analysis"},
    {"id": "אנליזה נומרית", "prereqs": ["מבוא מורחב למדמ״ח", "חדוא 2", "לינארית 2"], "coreqs": [], "category": "applications"},
    {"id": "חישוב מדעי", "prereqs": ["אנליזה נומרית", "מרוכבות 1"], "coreqs": [], "category": "applications"},
    {"id": "הילברט", "prereqs": ["חדוא 2", "לינארית 2", "מבוא מורחב למדמ״ח"], "coreqs": [], "category": "analysis"},
    {"id": "טופולוגיה", "prereqs": ["חדוא 2", "מבוא לתורת הקבוצות"], "coreqs": [], "category": "topology"},
    {"id": "אנליזה פונקציונלית", "prereqs": ["הילברט"], "coreqs": ["ממשיות", "טופולוגיה"], "category": "analysis"},
    {"id": "לוגיקה", "prereqs": ["מבוא לתורת הקבוצות"], "coreqs": [], "category": "foundations"},
    {"id": "אנליזה על יריעות", "prereqs": ["גאומטריה דיפרנציאלית", "טופולוגיה"], "coreqs": [], "category": "geometry"},
    {"id": "חדוא 4", "prereqs": ["חדוא 3"], "coreqs": []},
    {"id": "מד״ח", "prereqs": ["ODE"], "coreqs": ["חדוא 4","מרוכבות 1"]},
    {"id": "מרוכבות 2", "prereqs": ["מרוכבות 1"], "coreqs": []},
    {"id": "אנליזה הרמונית", "prereqs": ["ממשיות","מרוכבות 1"], "coreqs": []},
    {"id": "הסתברות למתמטיקאים", "prereqs": ["ממשיות","מבוא להסתברות"], "coreqs": []},
    {"id": "משחקי בורל", "prereqs": ["הסתברות למתמטיקאים"], "coreqs": []},
    {"id": "מבוא לתהליכים מקריים", "prereqs": ["הסתברות למתמטיקאים"], "coreqs": []},
    {"id": "תורת הקבוצות", "prereqs": ["לינארית 1", "חדוא 1"], "coreqs": []}
];

// Prepare the elements for Cytoscape
const elements = {
    nodes: [],
    edges: []
};

// Add nodes
courses.forEach(course => {
    elements.nodes.push({
        data: {
            id: course.id,
            label: course.id
        }
    });
});

// Add edges
courses.forEach(course => {
    // Add prerequisite edges
    course.prereqs.forEach(prereq => {
        elements.edges.push({
            data: {
                source: prereq,
                target: course.id,
                type: 'prereq'
            }
        });
    });

    // Add corequisite edges
    course.coreqs.forEach(coreq => {
        elements.edges.push({
            data: {
                source: coreq,
                target: course.id,
                type: 'coreq'
            }
        });
    });
});

// Initialize Cytoscape
const cy = cytoscape({
    container: document.getElementById('cy'),
    elements: elements,
    style: [
        {
            selector: 'node',
            style: {
                'label': 'data(label)',
                'text-valign': 'center',
                'text-halign': 'center',
                'background-color': '#b3d9ff',
                'text-wrap': 'wrap',
                'text-max-width': '100px',
                'font-size': '14px',
                'width': '120px',
                'height': '50px',
                'padding': '10px',
                'shape': 'roundrectangle',
                'border-width': '1px',
                'border-color': '#666',
                'transition-property': 'background-color, border-width, border-color, opacity, shadow-blur, shadow-color, font-size',
                'transition-duration': '0.2s',
                'shadow-blur': '0',
                'shadow-color': '#000',
                'shadow-opacity': 0
            }
        },
        {
            selector: 'edge[type="prereq"]',
            style: {
                'width': 2,
                'line-color': '#000',
                'target-arrow-color': '#000',
                'target-arrow-shape': 'triangle',
                'curve-style': 'bezier',
                'transition-property': 'opacity, line-color',
                'transition-duration': '0.2s'
            }
        },
        {
            selector: 'edge[type="coreq"]',
            style: {
                'width': 2,
                'line-color': '#666',
                'target-arrow-color': '#666',
                'target-arrow-shape': 'triangle',
                'line-style': 'dashed',
                'curve-style': 'bezier',
                'transition-property': 'opacity, line-color',
                'transition-duration': '0.2s'
            }
        }
    ],
    layout: {
        name: 'dagre',
        rankDir: 'TB',
        padding: 30,
        spacingFactor: 1.2,
        animate: false,
        rankSep: 60,
        nodeSep: 30,
        ranker: 'network-simplex'
    }
});

// Add zoom controls and limits
const minZoom = 0.3; // Minimum zoom level
const maxZoom = 2.5; // Maximum zoom level
const baseFontSize = 14; // Base font size in pixels

// Cache DOM elements and frequently used collections
const searchInput = document.getElementById('search');
let cachedNodes = cy.nodes();
let cachedEdges = cy.edges();

// Optimize the debounced search function further
const performSearch = _.debounce((searchTerm) => {
    if (!searchTerm) {
        resetToOriginalLayout();
        return;
    }

    searchTerm = searchTerm.toLowerCase();
    
    cy.startBatch();

    // Reset all nodes to default style first
    cachedNodes.style({
        'background-color': '#b3d9ff',
        'border-width': '1px',
        'border-color': '#666',
        'shadow-blur': '0',
        'shadow-opacity': 0
    });

    // Hide all elements
    cy.elements().addClass('hidden');

    // Find matching nodes
    const matchingNodes = cachedNodes.filter(node => 
        node.data('label').toLowerCase().includes(searchTerm)
    );

    if (matchingNodes.length > 0) {
        // Get only directly connected elements
        const directNeighbors = matchingNodes.neighborhood();
        const relevantElements = matchingNodes.union(directNeighbors);
        const relevantEdges = relevantElements.edgesWith(relevantElements);

        // Show relevant elements
        relevantElements.removeClass('hidden');
        relevantEdges.removeClass('hidden');

        // Ensure styles are applied in the correct order
        // First style the neighbors
        directNeighbors.nodes().forEach(node => {
            if (!matchingNodes.contains(node)) {  // Only style if not a matching node
                node.style({
                    'background-color': '#E3F2FD',
                    'border-width': '1.5px',
                    'border-color': '#1976D2',
                    'shadow-blur': '0',
                    'shadow-opacity': 0
                });
            }
        });

        // Then style the matching nodes (to ensure they take precedence)
        matchingNodes.forEach(node => {
            node.style({
                'background-color': '#4CAF50',
                'border-width': '2px',
                'border-color': '#2E7D32',
                'shadow-blur': '10',
                'shadow-opacity': 0.5
            });
        });

        // Optimize layout settings
        const searchLayout = {
            name: 'dagre',
            fit: true,
            padding: 35,
            animate: false,
            rankDir: 'TB',
            rankSep: 70,
            nodeSep: 50,
            ranker: 'network-simplex',
            spacingFactor: 1.0,
            edgeSep: 30,
            minLen: function(edge) {
                return edge.data('type') === 'prereq' ? 1.5 : 1;
            }
        };

        // Run layout only on visible elements
        const layout = relevantElements.layout(searchLayout);
        layout.run();

        // Ensure styles are maintained after layout
        requestAnimationFrame(() => {
            // Reapply styles after layout to ensure consistency
            directNeighbors.nodes().forEach(node => {
                if (!matchingNodes.contains(node)) {
                    node.style({
                        'background-color': '#E3F2FD',
                        'border-width': '1.5px',
                        'border-color': '#1976D2'
                    });
                }
            });

            matchingNodes.forEach(node => {
                node.style({
                    'background-color': '#4CAF50',
                    'border-width': '2px',
                    'border-color': '#2E7D32',
                    'shadow-blur': '10',
                    'shadow-opacity': 0.5
                });
            });
        });

        cy.animation({
            fit: {
                padding: 35
            },
            duration: 200,
            easing: 'ease-out'
        }).play().promise('completed').then(updateTextSize);
    } else {
        resetToOriginalLayout();
    }

    cy.endBatch();
}, 100);

// Add CSS classes for better performance
cy.style()
    .selector('.hidden')
    .style({
        'display': 'none',
        'opacity': 0
    })
    .update();

// Optimize resetToOriginalLayout function
function resetToOriginalLayout() {
    cy.startBatch();
    
    // Remove all custom styles efficiently
    cy.elements().removeClass('hidden');
    
    // Reset node styles in bulk
    cachedNodes.style({
        'background-color': '#b3d9ff',
        'border-width': '1px',
        'border-color': '#666',
        'shadow-blur': '0',
        'shadow-opacity': 0
    });

    // Reset edge styles efficiently
    cachedEdges.forEach(edge => {
        const type = edge.data('type');
        edge.style({
            'line-color': type === 'prereq' ? '#000' : '#666',
            'target-arrow-color': type === 'prereq' ? '#000' : '#666'
        });
    });

    // Optimize layout reset
    const layout = cy.layout({
        name: 'dagre',
        rankDir: 'TB',
        padding: 30,
        spacingFactor: 1.2,
        animate: false,
        rankSep: 60,
        nodeSep: 30,
        ranker: 'network-simplex'
    });

    layout.run();
    
    // Single animation for fitting
    cy.animation({
        fit: {
            padding: 30
        },
        duration: 150,
        easing: 'ease-out'
    }).play().promise('completed').then(() => {
        requestAnimationFrame(() => {
            cachedNodes.forEach(node => {
                originalPositions[node.id()] = node.position();
            });
            updateTextSize();
        });
    });

    cy.endBatch();
}

// Optimize text size updates
const updateTextSize = _.throttle(() => {
    const currentZoom = cy.zoom();
    const newFontSize = Math.max(
        Math.min(
            baseFontSize / currentZoom,
            baseFontSize * 1.5
        ),
        baseFontSize * 0.7
    );
    
    const newMaxWidth = Math.max(100 / currentZoom, 80);
    
    cy.startBatch();
    cachedNodes.style({
        'font-size': `${newFontSize}px`,
        'text-max-width': `${newMaxWidth}px`
    });
    cy.endBatch();
}, 50);

// Optimize event listeners
searchInput.addEventListener('input', (e) => {
    const searchTerm = e.target.value.trim();
    if (!searchTerm) {
        resetView(); // Use the same reset functionality when search is cleared
    } else {
        performSearch(searchTerm);
    }
}, { passive: true });

// Cache nodes and edges after initial layout
cy.on('layoutstop', () => {
    cachedNodes = cy.nodes();
    cachedEdges = cy.edges();
    updateTextSize();
});

// Modified zoom functions to include text size update
function zoomOut() {
    const currentZoom = cy.zoom();
    const newZoom = Math.max(currentZoom * 0.8, minZoom);
    cy.zoom({
        level: newZoom,
        renderedPosition: { x: cy.width() / 2, y: cy.height() / 2 }
    });
    updateTextSize();
}

function zoomIn() {
    const currentZoom = cy.zoom();
    const newZoom = Math.min(currentZoom * 1.2, maxZoom);
    cy.zoom({
        level: newZoom,
        renderedPosition: { x: cy.width() / 2, y: cy.height() / 2 }
    });
    updateTextSize();
}

// Add zoom event handler for continuous text size updates
cy.on('zoom', function() {
    updateTextSize();
});

function resetView() {
    // Clear search input and trigger search update
    if (searchInput) {
        searchInput.value = '';
        performSearch(''); // Explicitly trigger search with empty string
    }
    
    // Reset to original layout
    resetToOriginalLayout();
}

cy.minZoom(minZoom);
cy.maxZoom(maxZoom);

// Node tap and double click handlers
cy.on('tap', 'node', function(evt){
    const node = evt.target;
    cy.fit(node, 50);
});

cy.on('dblclick', function(evt){
    if(evt.target === cy){
        resetView();
    }
});

// Initial fit
cy.fit(50);

// Store the original layout for restoration
let originalPositions = {};

// After initial layout, store positions
cy.on('layoutstop', () => {
    cy.nodes().forEach(node => {
        originalPositions[node.id()] = { x: node.position('x'), y: node.position('y') };
    });
});

// Store the original layout configuration
const originalLayout = {
    name: 'dagre',
    rankDir: 'TB',
    padding: 50,
    spacingFactor: 1.5,
    animate: true,
    animationDuration: 500,
    rankSep: 100,
    nodeSep: 50,
    ranker: 'tight-tree'
};

// Modify the zoom handler
cy.on('zoom', function(evt) {
    const currentZoom = cy.zoom();
    if (currentZoom <= minZoom * 1.1) { // Restore when close to min zoom
        resetToOriginalLayout();
    }
    updateTextSize();
});

// Initial text size update
updateTextSize(); 