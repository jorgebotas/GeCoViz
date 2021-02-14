"use strict";
var __awaiter = (this && this.__awaiter) || function (thisArg, _arguments, P, generator) {
    function adopt(value) { return value instanceof P ? value : new P(function (resolve) { resolve(value); }); }
    return new (P || (P = Promise))(function (resolve, reject) {
        function fulfilled(value) { try { step(generator.next(value)); } catch (e) { reject(e); } }
        function rejected(value) { try { step(generator["throw"](value)); } catch (e) { reject(e); } }
        function step(result) { result.done ? resolve(result.value) : adopt(result.value).then(fulfilled, rejected); }
        step((generator = generator.apply(thisArg, _arguments || [])).next());
    });
};
var __generator = (this && this.__generator) || function (thisArg, body) {
    var _ = { label: 0, sent: function() { if (t[0] & 1) throw t[1]; return t[1]; }, trys: [], ops: [] }, f, y, t, g;
    return g = { next: verb(0), "throw": verb(1), "return": verb(2) }, typeof Symbol === "function" && (g[Symbol.iterator] = function() { return this; }), g;
    function verb(n) { return function (v) { return step([n, v]); }; }
    function step(op) {
        if (f) throw new TypeError("Generator is already executing.");
        while (_) try {
            if (f = 1, y && (t = op[0] & 2 ? y["return"] : op[0] ? y["throw"] || ((t = y["return"]) && t.call(y), 0) : y.next) && !(t = t.call(y, op[1])).done) return t;
            if (y = 0, t) op = [op[0] & 2, t.value];
            switch (op[0]) {
                case 0: case 1: t = op; break;
                case 4: _.label++; return { value: op[1], done: false };
                case 5: _.label++; y = op[1]; op = [0]; continue;
                case 7: op = _.ops.pop(); _.trys.pop(); continue;
                default:
                    if (!(t = _.trys, t = t.length > 0 && t[t.length - 1]) && (op[0] === 6 || op[0] === 2)) { _ = 0; continue; }
                    if (op[0] === 3 && (!t || (op[1] > t[0] && op[1] < t[3]))) { _.label = op[1]; break; }
                    if (op[0] === 6 && _.label < t[1]) { _.label = t[1]; t = op; break; }
                    if (t && _.label < t[2]) { _.label = t[2]; _.ops.push(op); break; }
                    if (t[2]) _.ops.pop();
                    _.trys.pop(); continue;
            }
            op = body.call(thisArg, _);
        } catch (e) { op = [6, e]; y = 0; } finally { f = t = 0; }
        if (op[0] & 5) throw op[1]; return { value: op[0] ? op[1] : void 0, done: true };
    }
};
exports.__esModule = true;
var d3_1 = require("d3");
function get_context(query, isCluster, cutoff, url) {
    return __awaiter(this, void 0, void 0, function () {
        var url_root, data, url_path;
        return __generator(this, function (_a) {
            switch (_a.label) {
                case 0:
                    url_root = "/gmgfam/getcontext/";
                    if (isCluster) {
                        url_root += "cluster/";
                    }
                    else {
                        url_root += "unigene/";
                    }
                    if (url) {
                        url_path = url;
                    }
                    else {
                        url_path = url_root + query + "/" + cutoff;
                    }
                    $('.loader').show();
                    return [4 /*yield*/, $.ajax({
                            url: url_path,
                            complete: function () {
                                $('.loader').hide();
                            },
                            error: function () {
                                $('.loader').hide();
                                alert("Incorrect cluster identifier");
                            }
                        })];
                case 1:
                    data = _a.sent();
                    return [2 /*return*/, data];
            }
        });
    });
}
function get_newick(url) {
    return __awaiter(this, void 0, void 0, function () {
        var newick;
        return __generator(this, function (_a) {
            switch (_a.label) {
                case 0: return [4 /*yield*/, $.ajax({
                        url: url,
                        error: function () {
                            console.log("No tree found");
                        }
                    })];
                case 1:
                    newick = _a.sent();
                    return [2 /*return*/, newick];
            }
        });
    });
}
function launch_analysis(selector, query, isCluster, cutoff, nenv, queryList) {
    return __awaiter(this, void 0, void 0, function () {
        var colors_file, colors, data, url, newick, newick_url, _a;
        return __generator(this, function (_b) {
            switch (_b.label) {
                case 0:
                    colors_file = "/static/geco/txt/colors.txt";
                    return [4 /*yield*/, fetch(colors_file)
                            .then(function (response) { return response.text(); })
                            .then(function (hex) { return colors = eval(hex); })
                        // ASYNC CALL TO MONGO SERVER
                    ];
                case 1:
                    _b.sent();
                    if (!queryList) return [3 /*break*/, 3];
                    url = '/getcontext/list/' + queryList.join(",") + "/" + cutoff;
                    return [4 /*yield*/, get_context(query, isCluster, cutoff, url)];
                case 2:
                    data = _b.sent();
                    return [3 /*break*/, 5];
                case 3: return [4 /*yield*/, get_context(query, isCluster, cutoff)];
                case 4:
                    data = _b.sent();
                    _b.label = 5;
                case 5:
                    _b.trys.push([5, 7, , 8]);
                    newick_url = "/gmgfam/tree/" + query + "/";
                    return [4 /*yield*/, get_newick(newick_url)];
                case 6:
                    newick = _b.sent();
                    return [3 /*break*/, 8];
                case 7:
                    _a = _b.sent();
                    newick = undefined;
                    return [3 /*break*/, 8];
                case 8:
                    // Launch GeCo which is a window global from geco.js
                    window.launch_GeCo(selector, data, newick, nenv, colors);
                    return [2 /*return*/];
            }
        });
    });
}
function gmgfam_geco(selector) {
    return __awaiter(this, void 0, void 0, function () {
        var nenv, urlPath, context, query, cutoff, genelist;
        return __generator(this, function (_a) {
            switch (_a.label) {
                case 0:
                    nenv = 41;
                    urlPath = window.location.pathname.split("/");
                    context = urlPath[2];
                    query = urlPath[3];
                    cutoff = +urlPath[4];
                    if (!(context == "clustercontext")) return [3 /*break*/, 2];
                    return [4 /*yield*/, launch_analysis(selector, query, true, cutoff, nenv)];
                case 1:
                    _a.sent();
                    return [3 /*break*/, 6];
                case 2:
                    if (!(context == "unigenecontext")) return [3 /*break*/, 4];
                    return [4 /*yield*/, launch_analysis(selector, query, false, cutoff, nenv)];
                case 3:
                    _a.sent();
                    return [3 /*break*/, 6];
                case 4:
                    if (!(context == "listcontext")) return [3 /*break*/, 6];
                    genelist = String(urlPath[3]).split(",");
                    return [4 /*yield*/, launch_analysis(selector, undefined, false, cutoff, nenv, genelist)];
                case 5:
                    _a.sent();
                    _a.label = 6;
                case 6:
                    d3_1.select('ul.navbar')
                        .style('visibility', 'visible')
                        .style('opacity', 1);
                    return [2 /*return*/];
            }
        });
    });
}
gmgfam_geco("body");
