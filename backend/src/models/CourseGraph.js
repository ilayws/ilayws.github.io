const { Graph } = require('graphlib');
const courseData = require('../data/courses.json');

class CourseGraph {
    constructor() {
        this.graph = new Graph({ directed: true });
        this.courseDetails = new Map();
        this.initialize();
    }

    initialize() {
        courseData.forEach(course => {
            this.addCourse(course);
        });
    }

    addCourse(course) {
        const { id, prereqs, coreqs, semester, difficulty, category } = course;
        
        // Add node if it doesn't exist
        if (!this.graph.hasNode(id)) {
            this.graph.setNode(id);
            this.courseDetails.set(id, { semester, difficulty, category, coreqs });
        }

        // Add prerequisite edges
        prereqs.forEach(prereq => {
            if (!this.graph.hasNode(prereq)) {
                this.graph.setNode(prereq);
            }
            this.graph.setEdge(prereq, id, 'prereq');
        });

        // Add corequisite edges (bidirectional)
        coreqs.forEach(coreq => {
            if (!this.graph.hasNode(coreq)) {
                this.graph.setNode(coreq);
            }
            this.graph.setEdge(coreq, id, 'coreq');
            this.graph.setEdge(id, coreq, 'coreq');
        });
    }

    getAllPrerequisites(courseId, includeIndirect = true) {
        if (!this.graph.hasNode(courseId)) {
            return [];
        }

        const prereqs = new Set();
        const visited = new Set();

        const dfs = (node) => {
            if (visited.has(node)) return;
            visited.add(node);

            this.graph.inEdges(node)?.forEach(edge => {
                if (edge.name === 'prereq') {
                    prereqs.add(edge.v);
                    if (includeIndirect) {
                        dfs(edge.v);
                    }
                }
            });
        };

        dfs(courseId);
        return Array.from(prereqs);
    }

    getDirectPrerequisites(courseId) {
        return this.getAllPrerequisites(courseId, false);
    }

    getCorequisites(courseId) {
        if (!this.graph.hasNode(courseId)) {
            return [];
        }

        const coreqs = new Set();
        this.graph.outEdges(courseId)?.forEach(edge => {
            if (edge.name === 'coreq') {
                coreqs.add(edge.w);
            }
        });

        return Array.from(coreqs);
    }

    getDependentCourses(courseId) {
        if (!this.graph.hasNode(courseId)) {
            return [];
        }

        const dependents = new Set();
        this.graph.outEdges(courseId)?.forEach(edge => {
            if (edge.name === 'prereq') {
                dependents.add(edge.w);
            }
        });

        return Array.from(dependents);
    }

    getPrereqChainDepth(courseId) {
        if (!this.graph.hasNode(courseId)) {
            return 0;
        }

        const visited = new Set();
        const dfs = (node) => {
            if (visited.has(node)) return 0;
            visited.add(node);

            let maxDepth = 0;
            this.graph.inEdges(node)?.forEach(edge => {
                if (edge.name === 'prereq') {
                    maxDepth = Math.max(maxDepth, 1 + dfs(edge.v));
                }
            });

            return maxDepth;
        };

        return dfs(courseId);
    }

    getAvailableCourses(completedCourses) {
        const available = new Set();
        
        this.graph.nodes().forEach(courseId => {
            const prereqs = this.getAllPrerequisites(courseId);
            const coreqs = this.getCorequisites(courseId);
            
            // Check if all prerequisites are completed
            const prereqsMet = prereqs.every(prereq => completedCourses.includes(prereq));
            
            // Check if at least one corequisite is completed or in progress
            const coreqsMet = coreqs.length === 0 || coreqs.some(coreq => 
                completedCourses.includes(coreq) || available.has(coreq)
            );

            if (prereqsMet && coreqsMet && !completedCourses.includes(courseId)) {
                available.add(courseId);
            }
        });

        return Array.from(available);
    }

    findShortestPath(fromCourseId, toCourseId) {
        if (!this.graph.hasNode(fromCourseId) || !this.graph.hasNode(toCourseId)) {
            return null;
        }

        const visited = new Set();
        const queue = [[fromCourseId, [fromCourseId]]];
        
        while (queue.length > 0) {
            const [current, path] = queue.shift();
            
            if (current === toCourseId) {
                return path;
            }

            if (!visited.has(current)) {
                visited.add(current);
                
                this.graph.outEdges(current)?.forEach(edge => {
                    if (!visited.has(edge.w)) {
                        queue.push([edge.w, [...path, edge.w]]);
                    }
                });
            }
        }

        return null;
    }

    getCourseDetails(courseId) {
        return this.courseDetails.get(courseId) || null;
    }

    getAllCourses() {
        return this.graph.nodes().map(courseId => ({
            id: courseId,
            ...this.getCourseDetails(courseId),
            prereqs: this.getDirectPrerequisites(courseId),
            coreqs: this.getCorequisites(courseId)
        }));
    }

    filterCourses({ semester, difficulty, category, hasPrereqs, hasCoreqs }) {
        return this.getAllCourses().filter(course => {
            if (semester && course.semester !== semester) return false;
            if (difficulty && course.difficulty !== difficulty) return false;
            if (category && course.category !== category) return false;
            if (hasPrereqs !== undefined) {
                const hasPrerequisites = this.getDirectPrerequisites(course.id).length > 0;
                if (hasPrereqs !== hasPrerequisites) return false;
            }
            if (hasCoreqs !== undefined) {
                const hasCorequisites = this.getCorequisites(course.id).length > 0;
                if (hasCoreqs !== hasCorequisites) return false;
            }
            return true;
        });
    }
}

module.exports = CourseGraph; 